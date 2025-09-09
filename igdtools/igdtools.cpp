/* Tools for converting to and manipulating IGD files.
 *
 * Run "igdtools --help" for usage.
 */
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_set>
#include <vector>

#include <sys/stat.h>

#include "picovcf.hpp"
#include "third-party/args.hxx"
#include "third-party/json.hpp"

using namespace picovcf;

const std::string DIR_SEP = "/";

// Supported metadata fields for conversion from VCF
enum MetadataFields {
    MDF_ALL = 0,
    MDF_CHROM = 1,
    MDF_QUAL = 2,
    MDF_FILTER = 3,
    MDF_INFO = 4,
};

// Supported meta-data fields for exporting. If you export "INFO" then you get all of INFO, there is no way
// to reduce it (use bcftools or something if you want that.)
static const std::string METADATA_FIELDS[] = {
    "ALL",
    "CHROM",
    "QUAL",
    "FILTER",
    "INFO",
};
// I want to maintain support for C++11, and std::size is not available until C++17, so we use
// a C-style array here.
#define ARRAY_SIZE(a) (sizeof(a) / sizeof(a[0]))
static_assert(ARRAY_SIZE(METADATA_FIELDS) == (size_t)MetadataFields::MDF_INFO + 1, "Enum/string mismatch");

/**
 * Conceptually similar to the Python multiprocessing.Pool, except using threads instead of processes.
 * All of the work must be added to the queue (via addWork) _prior_ to calling doAllWork. This is essentially
 * just a way to batch a bunch of known work.
 */
template <typename T, typename O> class PooledJobs {
public:
    PooledJobs() = default;
    virtual ~PooledJobs() = default;

    PooledJobs(PooledJobs& rhs) = delete;
    PooledJobs(PooledJobs&& rhs) = delete;
    PooledJobs& operator=(PooledJobs& rhs) = delete;
    PooledJobs& operator=(PooledJobs&& rhs) = delete;

    void addWork(T item) {
        std::lock_guard<std::mutex> lock(m_queueMutex);
        m_queue.push_back(std::move(item));
    }

    void doAllWork(size_t jobs) {
        if (jobs == 1) {
            while (!m_queue.empty()) {
                T nextItem = std::move(m_queue.front());
                m_queue.pop_front();
                m_results.push_back(processItem(std::move(nextItem)));
            }
        } else {
            std::vector<std::thread> threadPool;
            for (size_t i = 0; i < jobs; i++) {
                threadPool.emplace_back(&PooledJobs::workerMethod, this);
            }
            for (size_t i = 0; i < jobs; i++) {
                threadPool[i].join();
            }
        }
    }

    const std::list<O>& getResults() { return m_results; }

protected:
    // The worker (a single thread) loops until the queue is empty.
    void workerMethod() {
        while (true) {
            T nextItem; // T needs a default constructor.
            {
                std::lock_guard<std::mutex> lock(m_queueMutex);
                if (m_queue.empty()) {
                    break;
                }
                nextItem = std::move(m_queue.front());
                m_queue.pop_front();
            }
            O output = processItem(std::move(nextItem));
            {
                std::lock_guard<std::mutex> lock(m_queueMutex);
                m_results.push_back(std::move(output));
            }
        }
    }

    virtual O processItem(T item) = 0;

private:
    std::mutex m_queueMutex;
    std::list<T> m_queue;
    std::list<O> m_results;
};

inline void emitAllele(VariantT alleleIndex, std::ostream& out) {
    if (alleleIndex == MISSING_VALUE) {
        out << "? ";
    } else {
        out << alleleIndex << " ";
    }
}

inline bool ends_with(std::string const& string1, std::string const& string2) {
    if (string1.length() < string2.length()) {
        return false;
    }
    return (string2 == string1.substr(string1.length() - string2.length()));
}

template <typename Out> inline void split(const std::string& s, char delim, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

inline std::vector<std::string> split(const std::string& s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return std::move(elems);
}

inline std::string upper(const std::string& s) {
    std::string result = s;
    for (size_t i = 0; i < s.size(); i++) {
        result[i] = std::toupper(result[i]);
    }
    return std::move(result);
}

inline std::string lower(const std::string& s) {
    std::string result = s;
    for (size_t i = 0; i < s.size(); i++) {
        result[i] = std::tolower(result[i]);
    }
    return std::move(result);
}

inline std::string removeExt(const std::string& pathname) {
    size_t pos = pathname.find_last_of('.');
    if (pos != std::string::npos) {
        return std::move(pathname.substr(0, pos));
    }
    return pathname;
}

void make_dir_or_throw(const std::string& directory) {
    std::stringstream ssErr;
    int rv = mkdir(directory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (rv == 0) {
        return;
    }
    switch (errno) {
    case EEXIST: ssErr << "Directory " << directory << " already exists; remove and try again"; break;
    case EACCES: ssErr << "Access denied when creating directory " << directory; break;
    default: ssErr << "Failed creating directory " << directory << ": " << strerror(errno);
    }
    throw std::runtime_error(ssErr.str());
}

std::vector<std::string> readLines(const std::string& filename) {
    std::ifstream inStream(filename);
    if (!inStream.good()) {
        std::stringstream ssErr;
        ssErr << "Could not read file " << filename;
        throw std::runtime_error(ssErr.str());
    }
    std::string line;
    std::vector<std::string> result;
    while (std::getline(inStream, line)) {
        result.push_back(std::move(line));
    }
    return std::move(result);
}

// Supported types for VCF structured metadata (i.e., INFO fields)
enum VcfInfoType {
    VCFIT_INTEGER = 0,
    VCFIT_FLOAT = 1,
    VCFIT_STRING = 2,
};

// Just holds information for writing metadata info to output files.
struct MetadataWriteInfo {
    const std::string metadataOutDirectory;
    std::ofstream variantIdStream;
    std::unordered_map<std::string, std::ofstream> outStreams;
    std::unordered_map<std::string, std::ofstream> infoOutStreams;
    std::unordered_map<std::string, VcfInfoType> infoTypes;
    std::array<bool, ARRAY_SIZE(METADATA_FIELDS)> outputField;

    MetadataWriteInfo(const std::string& prefix, const std::string& metadataFieldList)
        : metadataOutDirectory(prefix) {

        make_dir_or_throw(metadataOutDirectory);

        for (size_t i = 0; i < outputField.size(); i++) {
            outputField[i] = false;
        }
        auto fields = split(metadataFieldList, ',');
        for (const auto& f : fields) {
            std::string fieldName = upper(f);
            bool matched = false;
            for (size_t i = 0; i < ARRAY_SIZE(METADATA_FIELDS); i++) {
                if (METADATA_FIELDS[i] == fieldName) {
                    outputField.at(i) = true;
                    matched = true;
                }
            }
            if (!matched) {
                PICOVCF_THROW_ERROR(
                    ApiMisuse,
                    "Unexpected metadata field name: " << f << ". Expected one of: all, chrom, qual, filter, info");
            }
        }
    }
};

// Write a variant-at-a-time metadata information.
void writeNextVariant(const VCFFile& vcfFile, VCFVariantView& variant, void* context) {
    PICOVCF_RELEASE_ASSERT(context != nullptr);
    MetadataWriteInfo* info = static_cast<MetadataWriteInfo*>(context);

    VCFVariantInfo allInfo = variant.parseToVariantInfo(/*basicInfoOnly=*/false);

    // One-time initialization
    if (info->outStreams.empty() && info->infoOutStreams.empty()) {
        // We always emit a file that contains the variant IDs, so that the metadata can be used
        // even if the IGD file is filtered downstream.
        std::stringstream variantFilename;
        variantFilename << info->metadataOutDirectory << DIR_SEP << "variants.txt";
        info->variantIdStream = std::ofstream(variantFilename.str());
        nlohmann::json metaInfo;
        metaInfo["type"] = "ID";
        info->variantIdStream << "# " << metaInfo << std::endl;

        auto infoSuffix = [&](const std::string& key) { return lower(METADATA_FIELDS[MDF_INFO]) + "." + key; };

        auto makeOutfile = [&](const std::string& suffix) {
            std::stringstream filename;
            filename << info->metadataOutDirectory << DIR_SEP << suffix << ".txt";
            std::ofstream outStream(filename.str());
            if (!outStream) {
                std::stringstream ssErr;
                ssErr << "Cannot create output file " << filename.str();
                throw std::runtime_error(ssErr.str());
            }
            return std::move(outStream);
        };

        if (info->outputField[MDF_ALL] || info->outputField[MDF_CHROM]) {
            const auto key = METADATA_FIELDS[MDF_CHROM];
            nlohmann::json metaInfo;
            metaInfo["type"] = "CHROM";
            info->outStreams.emplace(key, std::move(makeOutfile("chrom")));
            info->outStreams.at(key) << "# " << metaInfo << std::endl;
        }
        if (info->outputField[MDF_ALL] || info->outputField[MDF_QUAL]) {
            const auto key = METADATA_FIELDS[MDF_QUAL];
            nlohmann::json metaInfo;
            metaInfo["type"] = "QUAL";
            info->outStreams.emplace(key, std::move(makeOutfile("qual")));
            info->outStreams.at(key) << "# " << metaInfo << std::endl;
        }
        if (info->outputField[MDF_ALL] || info->outputField[MDF_FILTER]) {
            const auto key = METADATA_FIELDS[MDF_FILTER];
            nlohmann::json metaInfo;
            metaInfo["type"] = "FILTER";
            info->outStreams.emplace(METADATA_FIELDS[MDF_FILTER], std::move(makeOutfile("filter")));
            info->outStreams.at(key) << "# " << metaInfo << std::endl;
        }
        if (info->outputField[MDF_ALL] || info->outputField[MDF_INFO]) {
            auto metaStrings = vcfFile.getAllMetaInfo(METADATA_FIELDS[MDF_INFO].c_str());
            for (auto s : metaStrings) {
                nlohmann::json metaInfo = picovcf_parse_structured_meta(s);
                metaInfo["type"] = "INFO";
                const auto suffix = infoSuffix(metaInfo["ID"]);
                const auto key = metaInfo["ID"];
                info->infoOutStreams.emplace(key, std::move(makeOutfile(suffix)));
                info->infoOutStreams.at(key) << "# " << metaInfo << std::endl;
                PICOVCF_ASSERT_OR_MALFORMED(metaInfo.contains("Type"), "INFO field is missing a \"Type\" key/value");
                VcfInfoType infoType = VCFIT_STRING;
                if (metaInfo["Type"] == "Integer") {
                    infoType = VCFIT_INTEGER;
                } else if (metaInfo["Type"] == "Float") {
                    infoType = VCFIT_FLOAT;
                } else if (metaInfo["Type"] != "String") {
                    PICOVCF_ASSERT_OR_MALFORMED(false, "Unknown metadata INFO type: " << metaInfo["Type"]);
                }
                info->infoTypes.emplace(key, infoType);
            }
        }
    }

    info->variantIdStream << variant.getID() << std::endl;

    if (info->outputField[MDF_ALL] || info->outputField[MDF_CHROM]) {
        info->outStreams.at(METADATA_FIELDS[MDF_CHROM]) << allInfo.chromosome << std::endl;
    }
    if (info->outputField[MDF_ALL] || info->outputField[MDF_QUAL]) {
        info->outStreams.at(METADATA_FIELDS[MDF_QUAL]) << allInfo.quality << std::endl;
    }
    if (info->outputField[MDF_ALL] || info->outputField[MDF_FILTER]) {
        info->outStreams.at(METADATA_FIELDS[MDF_FILTER]) << allInfo.filter << std::endl;
    }
    if (info->outputField[MDF_ALL] || info->outputField[MDF_INFO]) {
        for (const auto& keyStream : info->infoOutStreams) {
            auto infoIt = allInfo.information.find(keyStream.first);
            if (infoIt != allInfo.information.end()) {
                info->infoOutStreams.at(keyStream.first) << infoIt->second << std::endl;
            } else {
                switch (info->infoTypes.at(keyStream.first)) {
                case VCFIT_INTEGER: info->infoOutStreams.at(keyStream.first) << "0" << std::endl; break;
                case VCFIT_FLOAT: info->infoOutStreams.at(keyStream.first) << "NaN" << std::endl; break;
                case VCFIT_STRING: info->infoOutStreams.at(keyStream.first) << "." << std::endl; break;
                default: PICOVCF_RELEASE_ASSERT(false);
                }
            }
        }
    }
}

struct VCFExportMetaJob {
    picovcf::RangePair range;
    std::string vcfFilename;
    std::string metadataOutDirectory;
    std::string metadataFieldList;
};

class VCFExportWorker : public PooledJobs<VCFExportMetaJob, bool> {
public:
    // Given a VCF file, just emit the metadata from it, don't convert it to IGD.
    static void vcfExportMetadata(const picovcf::RangePair& range,
                                  const std::string& vcfFilename,
                                  const std::string& metadataOutDirectory,
                                  const std::string& metadataFieldList) {
        VCFFile vcf(vcfFilename);

        MetadataWriteInfo writeInfo(metadataOutDirectory, metadataFieldList);

        bool alreadyHave = false;
        if (range.first == INVALID_POSITION) {
            vcf.seekBeforeVariants();
        } else {
            alreadyHave = vcf.lowerBoundPosition(range.first);
        }
        while (alreadyHave || vcf.nextVariant()) {
            alreadyHave = false;
            VCFVariantView& variant = vcf.currentVariant();
            const auto position = variant.getPosition();
            if (position > range.second) {
                break;
            }
            for (size_t i = 0; i < variant.getAltAlleles().size(); i++) {
                writeNextVariant(vcf, variant, &writeInfo);
            }
        }
    }

protected:
    bool processItem(VCFExportMetaJob job) override {
        vcfExportMetadata(job.range, job.vcfFilename, job.metadataOutDirectory, job.metadataFieldList);
        return true;
    }
};

struct VCFConvertJob {
    size_t jobId;
    std::string vcfFilename;
    std::string igdFilename;
    std::string description;
    std::string exportList;
    bool emitVariantIds;
    bool forceUnphased;
    size_t forceToPloidy;
    bool dropUnphased;
    picovcf::RangePair range;
};

class VCFConvertWorker : public PooledJobs<VCFConvertJob, std::string> {
public:
    static std::string vcfConvert(const VCFConvertJob& job) {
        void (*variantCallback)(const VCFFile&, VCFVariantView&, void*) = nullptr;
        void* callbackContext = nullptr;

        std::unique_ptr<MetadataWriteInfo> writeInfo;
        if (!job.exportList.empty()) {
            std::string metadataOutDirectory = job.igdFilename;
            if (ends_with(metadataOutDirectory, ".igd")) {
                metadataOutDirectory = removeExt(metadataOutDirectory);
            }
            std::stringstream ssMDOut;
            ssMDOut << metadataOutDirectory << ".meta";
            if (job.jobId > 0) {
                ssMDOut << job.jobId;
            }
            std::cout << "Exporting metadata to directory: " << ssMDOut.str() << std::endl;

            writeInfo.reset(new MetadataWriteInfo(ssMDOut.str(), job.exportList));
            variantCallback = writeNextVariant;
            callbackContext = writeInfo.get();
        }
        std::stringstream outFilename;
        outFilename << job.igdFilename;
        if (job.jobId > 0) {
            outFilename << "." << job.jobId;
        }
        vcfToIGD(job.vcfFilename,
                 outFilename.str(),
                 job.description,
                 true,
                 /*emitIndividualIds=*/true,
                 job.emitVariantIds,
                 job.forceUnphased,
                 job.forceToPloidy,
                 job.dropUnphased,
                 variantCallback,
                 callbackContext,
                 job.range);
        return outFilename.str();
    }

protected:
    std::string processItem(VCFConvertJob job) override { return vcfConvert(job); }
};

std::vector<picovcf::RangePair> getRangesForJobs(const std::string& vcfFilename, size_t numJobs) {
    picovcf::VCFFile vcf(vcfFilename);
    picovcf::RangePair fileBPRange = vcf.getGenomeRange();
    const size_t span = fileBPRange.second - fileBPRange.first;
    const size_t perPart = (span + (numJobs - 1)) / numJobs;
    if (perPart < 100 || numJobs == 1) {
        return {fileBPRange};
    }
    if (!vcf.isUsingIndex()) {
        std::cerr << "WARNING: Parallel processing a VCF file without a tabix index is NOT recommended" << std::endl;
    }
    std::vector<picovcf::RangePair> result;
    size_t start = fileBPRange.first;
    while (start <= fileBPRange.second) {
        result.emplace_back(start, std::min(fileBPRange.second, start + perPart - 1));
        start += perPart;
    }
    return std::move(result);
}

int main(int argc, char* argv[]) {
    args::ArgumentParser parser("Process or create IGD files.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::Positional<std::string> infile(parser, "input_file", "The input file (.vcf, .vcf.gz, or .igd)");
    args::ValueFlag<std::string> outfile(parser, "output", "The output file to produce.", {'o', "out"});
    args::ValueFlag<std::string> outputDescription(
        parser, "outputDescription", "The description string to include in the IGD output.", {"description"});
    args::ValueFlagList<std::string> mergeWith(
        parser, "mergeWith", "Merge the input IGD file with all of these specified IGD files.", {"merge"});
    args::ValueFlag<std::string> range(
        parser, "range", "Restrict to the given base-pair range (inclusive).", {'r', "range"});
    args::ValueFlag<std::string> frange(parser,
                                        "frange",
                                        "Restrict to variants with frequency in "
                                        "the range (inclusive, exclusive).",
                                        {'f', "frange"});
    args::Flag info(parser, "info", "Display information from the IGD header.", {'i', "info"});
    args::Flag individuals(parser, "individuals", "Emit the mapping from individual index to ID.", {"individuals"});
    args::Flag variants(parser, "variants", "Emit the mapping from variant index to ID.", {"variants"});
    args::Flag stats(parser,
                     "stats",
                     "Emit some simple statistics about the distribution of "
                     "samples, sites, and variants.",
                     {'s', "stats"});
    args::Flag alleles(parser, "alleles", "Emit allele frequencies.", {'a', "alleles"});
    args::Flag noVariantIds(
        parser, "noVariantIds", "Do not emit IDs for variants in the resulting IGD file.", {"no-var-ids"});
    args::ValueFlag<size_t> trimSamples(parser, "trimSamples", "Trim samples to this many individuals.", {"trim"});
    args::ValueFlag<std::string> keepSamples(
        parser, "keepSamples", "Keep samples listed in the given text file.", {'S', "samples"});
    args::Flag forceUnphasedArg(
        parser, "forceUnphased", "Force output file to be unphased, regardless of input.", {"force-unphased"});
    args::Flag dropMultiSites(
        parser, "dropMultiSites", "Drop multi-allelic sites (more than one alternate allele).", {"drop-multi-sites"});
    args::Flag dropNonSNVs(
        parser, "dropNonSNVs", "Drop variants containing alleles that are not single nucleotides.", {"drop-non-snvs"});
    args::Flag dropNonSNVSites(parser,
                               "dropNonSNVSites",
                               "Drop sites containing alleles that are not single nucleotides.",
                               {"drop-non-snv-sites"});
    args::Flag dropUnphased(parser, "dropUnphased", "Drop sites containing unphased data.", {"drop-unphased"});
    args::ValueFlag<std::string> exportMetadata(
        parser,
        "exportMetadata",
        "Export the metadata from the given .vcf[.gz] file to the given filename. The output format is a "
        "text file in the format that can be loaded via numpy.loadtxt(), where the first line is a comment "
        "containing information about the metadata. This option takes a list of metadata to export, which "
        "can be: all, chrom, qual, filter, info",
        {'e', "export-metadata"});
    args::ValueFlag<size_t> forceToPloidy(parser,
                                          "forceToPloidy",
                                          "IGD files have a single ploidy, but you can use this to force all "
                                          "variants/individuals to have the same ploidy.",
                                          {'p', "force-ploidy"});
    args::ValueFlag<std::string> updateIndividualIds(
        parser,
        "update-indiv-ids",
        "Given a text file with an individual ID per line, set the IDs in the output IGD file to match them",
        {"update-indiv-ids"});
    args::ValueFlag<size_t> jobs(parser,
                                 "jobs",
                                 "How many jobs (threads) to use. Only VCF->IGD conversion supports this currently",
                                 {'j', "jobs"});
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help&) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    } catch (args::ValidationError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (!infile) {
        std::cout << parser;
        return 0;
    }

#define UNSUPPORTED_FOR_VCF(parameter, parameterName)                                                                  \
    do {                                                                                                               \
        if ((bool)(parameter)) {                                                                                       \
            std::cerr << "Parameter " << (parameterName) << " is not supported for VCF file conversion" << std::endl;  \
            return 1;                                                                                                  \
        }                                                                                                              \
    } while (0)

#define UNSUPPORTED_FOR_MERGE(parameter, parameterName)                                                                \
    do {                                                                                                               \
        if ((bool)(parameter)) {                                                                                       \
            std::cerr << "Parameter " << (parameterName) << " is not supported for merging" << std::endl;              \
            return 1;                                                                                                  \
        }                                                                                                              \
    } while (0)

#define ONLY_SUPPORTED_FOR_VCF(parameter, parameterName)                                                               \
    do {                                                                                                               \
        if ((bool)(parameter)) {                                                                                       \
            std::cerr << "Parameter " << (parameterName) << " is only supported for VCF file conversion" << std::endl; \
            return 1;                                                                                                  \
        }                                                                                                              \
    } while (0)

    std::string description = outputDescription ? *outputDescription : "";

    try {
        size_t bpStart = 0;
        size_t bpEnd = std::numeric_limits<size_t>::max();
        if (range) {
            const char* rangec = range->c_str();
            const char* rangecEnd = range->c_str() + range->size();
            char* endptr = nullptr;
            bpStart = std::strtoull(rangec, &endptr, 10);
            if (endptr >= rangecEnd || *endptr != '-') {
                std::cerr << "Malformed range: " << *range << " (must be \"lower-upper\")" << std::endl;
                return 1;
            }
            endptr++;
            bpEnd = std::strtoull(endptr, &endptr, 10);
            if (endptr != rangecEnd) {
                std::cerr << "Malformed range: " << *range << " (must be \"lower-upper\")" << std::endl;
                return 1;
            }
            std::cout << "Restricting to base-pair range " << bpStart << " - " << bpEnd << " (inclusive)" << std::endl;
        }

        const bool isVcf = ends_with(*infile, ".vcf") || ends_with(*infile, ".vcf.gz");
        if (isVcf) {
            constexpr size_t MAX_JOBS = 1024;
            const size_t numJobs = jobs ? *jobs : 1;
            if (numJobs < 1 || numJobs > MAX_JOBS) {
                PICOVCF_THROW_ERROR(picovcf::ApiMisuse, "Invalid number of jobs specified");
            }

            // TODO: user-input range should apply to this as well...
            std::vector<picovcf::RangePair> jobRanges = getRangesForJobs(*infile, numJobs);
            PICOVCF_RELEASE_ASSERT(!jobRanges.empty());

            // FIXME: apply range to metadata export as well.
            if (!outfile) {
                if (exportMetadata) {
                    std::string metadataOutDirectory = *infile;
                    if (ends_with(metadataOutDirectory, ".gz")) {
                        metadataOutDirectory = removeExt(metadataOutDirectory);
                    }
                    if (ends_with(metadataOutDirectory, ".vcf")) {
                        metadataOutDirectory = removeExt(metadataOutDirectory);
                    }
                    metadataOutDirectory = metadataOutDirectory + ".meta";
                    std::cout << "Exporting metadata to file(s) beginning with prefix " << metadataOutDirectory
                              << std::endl;
                    VCFExportWorker exportWorker;
                    if (jobRanges.size() == 1) {
                        exportWorker.vcfExportMetadata(
                            jobRanges.front(), *infile, metadataOutDirectory, *exportMetadata);
                    } else {
                        size_t jobId = 1;
                        for (const auto& range : jobRanges) {
                            std::stringstream metaOutJob;
                            metaOutJob << metadataOutDirectory << jobId++;
                            std::cerr << "ADD_WORK(" << range.first << ", " << metaOutJob.str() << ")\n";
                            exportWorker.addWork({range, *infile, metaOutJob.str(), *exportMetadata});
                        }
                        exportWorker.doAllWork(numJobs);
                    }
                    return 0;
                }
                std::cerr << "VCF input is only supported for conversion to IGD, and/or to export meta-data. "
                          << "Use --out to convert to IGD, and --export-metadata to export the metadata." << std::endl;
                return 1;
            }
            UNSUPPORTED_FOR_VCF(dropMultiSites, "--drop-multi-sites");
            UNSUPPORTED_FOR_VCF(dropNonSNVs, "--drop-non-snvs");
            UNSUPPORTED_FOR_VCF(dropNonSNVSites, "--drop-non-snv-sites");
            UNSUPPORTED_FOR_VCF(keepSamples, "--samples");
            UNSUPPORTED_FOR_VCF(mergeWith, "--merge");

            // TODO: move all of this into a function, and then:
            // 1. Check if VCF is indexed/indexable -- vcf.gz and have HTS support
            //      --> if not, reject the -j parameter
            // 2. Come up with list of ranges [lower, upper) -- is this fast with an indexed VCF?
            // 3. Create a tmp directory for each part
            // 4. Convert that part to IGD (threaded)
            // 5. Call IGD merge (C++ impl needed)

            const bool emitVariantIds = !noVariantIds;

            VCFConvertWorker convertWorker;
            if (jobRanges.size() == 1) {
                convertWorker.vcfConvert({
                    /*jobId=*/0,
                    /*vcfFilename=*/*infile,
                    /*igdFilename=*/*outfile,
                    /*description=*/description,
                    /*exportList=*/exportMetadata ? *exportMetadata : "",
                    /*emitVariantIds=*/emitVariantIds,
                    /*forceUnphased=*/forceUnphasedArg,
                    /*forceToPloidy=*/forceToPloidy ? *forceToPloidy : 0,
                    /*dropUnphased=*/dropUnphased,
                    /*range=*/jobRanges.front(),
                });
            } else {
                size_t jobId = 1;
                for (const auto& range : jobRanges) {
                    VCFConvertJob job = {/*jobId=*/jobId++,
                                         /*vcfFilename=*/*infile,
                                         /*igdFilename=*/*outfile,
                                         /*description=*/description,
                                         /*exportList=*/exportMetadata ? *exportMetadata : "",
                                         /*emitVariantIds=*/emitVariantIds,
                                         /*forceUnphased=*/forceUnphasedArg,
                                         /*forceToPloidy=*/forceToPloidy ? *forceToPloidy : 0,
                                         /*dropUnphased=*/dropUnphased,
                                         /*range=*/range};
                    convertWorker.addWork(job);
                }
                convertWorker.doAllWork(numJobs);
                std::vector<std::string> mergeInputs;
                for (const auto& igdName : convertWorker.getResults()) {
                    mergeInputs.emplace_back(igdName);
                }
                std::ofstream outStream(*outfile, std::ios::binary);
                mergeIGDs(outStream, mergeInputs);
            }
            return 0;
        } else {
            ONLY_SUPPORTED_FOR_VCF(forceToPloidy, "--force-ploidy");
            ONLY_SUPPORTED_FOR_VCF(dropUnphased, "--drop-unphased");
            ONLY_SUPPORTED_FOR_VCF(exportMetadata, "--export-metadata");
            ONLY_SUPPORTED_FOR_VCF(jobs, "--jobs");
        }

        if (mergeWith) {
            if (!outfile) {
                std::cerr << "ERROR: --out (-o) required for --merge" << std::endl;
                return 2;
            }
            UNSUPPORTED_FOR_MERGE(range, "--range");
            UNSUPPORTED_FOR_MERGE(frange, "--frange");
            UNSUPPORTED_FOR_MERGE(info, "--info");
            UNSUPPORTED_FOR_MERGE(individuals, "--individuals");
            UNSUPPORTED_FOR_MERGE(variants, "--variants");
            UNSUPPORTED_FOR_MERGE(stats, "--stats");
            UNSUPPORTED_FOR_MERGE(alleles, "--alleles");
            UNSUPPORTED_FOR_MERGE(trimSamples, "--trim");
            UNSUPPORTED_FOR_MERGE(keepSamples, "--samples");
            UNSUPPORTED_FOR_MERGE(forceUnphasedArg, "--force-unphased");
            UNSUPPORTED_FOR_MERGE(dropMultiSites, "--drop-multi-sites");
            UNSUPPORTED_FOR_MERGE(dropNonSNVs, "--drop-non-snvs");
            UNSUPPORTED_FOR_MERGE(dropNonSNVSites, "--drop-non-snv-sites");
            UNSUPPORTED_FOR_MERGE(updateIndividualIds, "--update-individuals");
            UNSUPPORTED_FOR_MERGE(noVariantIds, "--no-var-ids");

            std::vector<std::string> inputFilenames = {*infile};
            for (const auto& fn : *mergeWith) {
                inputFilenames.emplace_back(fn);
            }
            std::ofstream outStream(*outfile, std::ios::binary);
            mergeIGDs(outStream, inputFilenames, outputDescription ? *outputDescription : "");
            return 0;
        }

        // Not a VCF, then assume it is IGD and load the header.
        IGDData igd(*infile);

        double fLower = 0.0;
        double fUpper = 1.1;
        if (frange) {
            const char* frangec = frange->c_str();
            const char* frangecEnd = frange->c_str() + frange->size();
            char* endptr = nullptr;
            fLower = std::strtod(frangec, &endptr);
            if (endptr >= frangecEnd || *endptr != '-') {
                std::cerr << "Malformed frange: " << *frange << " (must be \"lower-upper\")" << std::endl;
                return 1;
            }
            endptr++;
            fUpper = std::strtod(endptr, &endptr);
            if (endptr != frangecEnd) {
                std::cerr << "Malformed frange: " << *frange << " (must be \"lower-upper\")" << std::endl;
                return 1;
            }
            std::cout << "Restricting to allele frequencies between [" << fLower << ", " << fUpper << ")" << std::endl;
        }

        const uint64_t ploidy = igd.getPloidy();
        const size_t numSamples = igd.numSamples();
        size_t effectiveSampleCt = numSamples;

        if (info) {
            std::cout << "Header information for " << *infile << std::endl;
            if (igd.isPhased()) {
                std::cout << "  Variants: " << igd.numVariants() << std::endl;
            } else {
                std::cout << "  Variant rows (includes multiple rows for different num_copies): " << igd.numVariants()
                          << std::endl;
                std::cout << "  Variants: between " << igd.numVariants() / igd.getPloidy() << " and "
                          << igd.numVariants() << std::endl;
            }
            std::cout << "  Individuals: " << igd.numIndividuals() << std::endl;
            std::cout << "  Ploidy: " << igd.getPloidy() << std::endl;
            std::cout << "  Phased?: " << (igd.isPhased() ? "true" : "false") << std::endl;
            std::cout << "  Source: " << igd.getSource() << std::endl;
            std::cout << "  Genome range: " << igd.getPosition(0) << "-" << igd.getPosition(igd.numVariants() - 1)
                      << std::endl;
            std::cout << "  Has individual IDs? " << (igd.getIndividualIds().empty() ? "No" : "Yes") << std::endl;
            std::cout << "  Has variant IDs? " << (igd.getVariantIds().empty() ? "No" : "Yes") << std::endl;
        }
        if (individuals) {
            std::vector<std::string> individualIds = igd.getIndividualIds();
            if (individualIds.empty()) {
                std::cout << "No individual IDs in this IGD file" << std::endl;
            }
            for (size_t i = 0; i < individualIds.size(); i++) {
                std::cout << i << ": " << individualIds[i] << std::endl;
            }
        }
        if (variants) {
            std::vector<std::string> variantIds = igd.getVariantIds();
            if (variantIds.empty()) {
                std::cout << "No variant IDs in this IGD file" << std::endl;
            }
            for (size_t i = 0; i < variantIds.size(); i++) {
                std::cout << i << ": " << variantIds[i] << std::endl;
            }
        }

#define CONDITION_PRINT(condition, message)                                                                            \
    do {                                                                                                               \
        if (condition) {                                                                                               \
            std::cout << message;                                                                                      \
        }                                                                                                              \
    } while (0)

        const bool forcingUnphased = forceUnphasedArg && igd.isPhased();
        if (forcingUnphased && ploidy > std::numeric_limits<uint8_t>::max()) {
            std::cerr << "Cannot store data with ploidy " << ploidy << " as unphased" << std::endl;
            return 2;
        }

        // Drop any individuals that are not kept, and map them to the order in the filter list.
        std::vector<size_t> mapIndividual(igd.numIndividuals());
        std::iota(mapIndividual.begin(), mapIndividual.end(), 0);
        if (keepSamples) {
            if (trimSamples) {
                throw std::runtime_error("Cannot specify both --samples and --trim");
            }
            const auto& keepIndivs = readLines(*keepSamples);
            std::unordered_map<std::string, size_t> keepMap;
            size_t i = 0;
            for (const auto& indiv : keepIndivs) {
                if (!indiv.empty()) {
                    keepMap.emplace(indiv, i++);
                }
            }
            if (keepMap.empty()) {
                throw std::runtime_error("Provided sample list is empty");
            }
            size_t discarded = 0;
            const auto& allIndivs = igd.getIndividualIds();
            for (size_t i = 0; i < allIndivs.size(); i++) {
                const auto& findIt = keepMap.find(allIndivs[i]);
                if (findIt == keepMap.end()) {
                    mapIndividual.at(i) = std::numeric_limits<size_t>::max();
                    discarded++;
                } else {
                    mapIndividual.at(i) = findIt->second;
                }
            }
            const size_t keptIndivs = igd.numIndividuals() - discarded;
            if (keptIndivs != keepMap.size()) {
                std::cerr << "ERROR: Only kept " << keptIndivs << ", but sample list had " << keepMap.size()
                          << std::endl;
                return 2;
            }
            effectiveSampleCt = igd.isPhased() ? (ploidy * (keptIndivs)) : keptIndivs;
        } else if (trimSamples) {
            if (*trimSamples >= igd.numIndividuals()) {
                std::cerr << "--trim value is larger than or equal to the number of individuals" << std::endl;
                return 2;
            }
            for (size_t i = *trimSamples; i < igd.numIndividuals(); i++) {
                mapIndividual.at(i) = std::numeric_limits<size_t>::max();
            }
            effectiveSampleCt = igd.isPhased() ? (ploidy * (*trimSamples)) : *trimSamples;
        }

        std::shared_ptr<std::ofstream> igdOutfile;
        std::shared_ptr<IGDWriter> writer;
        const uint64_t numIndividuals = igd.isPhased() ? (effectiveSampleCt / ploidy) : effectiveSampleCt;
        if (outfile) {
            igdOutfile = std::make_shared<std::ofstream>(*outfile, std::ios::binary);
            writer = std::make_shared<IGDWriter>(ploidy, numIndividuals, !forcingUnphased && igd.isPhased());
            writer->writeHeader(*igdOutfile, *infile, igd.getDescription());
        }

        const bool iterateSamples = (bool)stats || (bool)alleles || (bool)outfile;
        if (iterateSamples) {
            static constexpr char SEP = '\t';
            std::stringstream copySS;
            if (!igd.isPhased()) {
                copySS << "COPIES" << SEP;
            }
            CONDITION_PRINT(alleles,
                            "POSITION" << SEP << copySS.str() << "REF" << SEP << "ALT" << SEP << "ALT COUNT" << SEP
                                       << "TOTAL" << std::endl);

            size_t sites = 0;
            bool _ignore = false;
            size_t variants = 0;
            size_t missingRows = 0;
            size_t missingAlleles = 0;
            std::vector<size_t> samplesPerVariant;
            size_t sampleRefsTotal = 0;
            std::vector<size_t> sampleToMuts(effectiveSampleCt); // Counts muts per sample

            std::vector<bool> skipVariant(igd.numVariants(), false);
            std::unordered_set<size_t> skipPosition;
            size_t skippedVars = 0;
            // Drop multi-allelic sites if requested.
            if (dropMultiSites) {
                size_t lastPosition = std::numeric_limits<size_t>::max();
                for (size_t i = 0; i < igd.numVariants(); i++) {
                    uint8_t numCopies = 0;
                    bool isMissing = false;
                    auto pos = igd.getPosition(i, isMissing, numCopies);
                    // Ignore missing data rows, they aren't relevant to this calculation.
                    if (isMissing) {
                        continue;
                    }
                    if (pos == lastPosition) {
                        skipPosition.emplace(pos);
                    }
                    lastPosition = pos;
                }
            }

            // Drop non-SNV sites/variants if requested.
            if (dropNonSNVSites || dropNonSNVs) {
                for (size_t i = 0; i < igd.numVariants(); i++) {
                    uint8_t numCopies = 0;
                    bool isMissing = false;
                    const size_t pos = igd.getPosition(i, isMissing, numCopies);
                    const bool nonSNV = (igd.getRefAllele(i).size() >= 2) || (igd.getAltAllele(i).size() >= 2);
                    if (nonSNV) {
                        // Entire site can be dropped, or just the variant, depending on user option chosen.
                        if (dropNonSNVSites) {
                            skipPosition.emplace(pos);
                        }
                        if (dropNonSNVs) {
                            skipVariant[i] = true;
                            skippedVars++;
                        }
                    }
                }
            }

            const size_t iDivisor = igd.isPhased() ? ploidy : 1;
            auto filterSamples = [&](IGDSampleList& samples) {
                size_t k = 0;
                for (size_t j = 0; j < samples.size(); j++) {
                    const size_t indiv = samples[j] / iDivisor;
                    size_t hap = 0;
                    if (igd.isPhased()) {
                        hap = samples[j] % iDivisor;
                    }
                    const size_t mappedIndex = mapIndividual.at(indiv);
                    if (mappedIndex != std::numeric_limits<size_t>::max()) {
                        samples[k++] = (mappedIndex * iDivisor) + hap;
                    }
                }
                if (k != samples.size()) {
                    samples.resize(k);
                }
            };

            if (!skipPosition.empty()) {
                std::cerr << "Skipping " << skipPosition.size() << " sites due to filters" << std::endl;
            }
            if (skippedVars > 0) {
                std::cerr << "Skipping " << skippedVars << " variants due to filters" << std::endl;
            }

            // Number of variants, counted by by num_copies
            std::vector<size_t> variantsByCopy(igd.getPloidy() + 1, 0);

            auto vids = igd.getVariantIds();
            std::vector<std::string> newVariantIds;
            size_t lastPosition = std::numeric_limits<size_t>::max();
            for (size_t i = 0; i < igd.numVariants(); i++) {
                uint8_t numCopies = 0;
                bool isMissing = false;
                auto pos = igd.getPosition(i, isMissing, numCopies);
                auto skipIt = skipPosition.find(pos);
                if (skipIt != skipPosition.end() || skipVariant[i]) {
                    continue;
                }
                if (pos >= bpStart && pos <= bpEnd) {
                    auto sampleList = igd.getSamplesWithAlt(i);
                    filterSamples(sampleList);
                    const auto sampleCt = sampleList.size();
                    assert(sampleCt <= effectiveSampleCt);
                    const double freq = (double)sampleCt / (double)effectiveSampleCt;
                    if (freq < fLower || freq >= fUpper) {
                        continue;
                    }
                    const auto& ref = igd.getRefAllele(i);
                    const auto& alt = igd.getAltAllele(i);
                    if (alleles) {
                        std::cout << pos << SEP;
                        if (!igd.isPhased()) {
                            std::cout << (uint64_t)numCopies << SEP;
                        }
                        const std::string altOut = isMissing ? "." : alt;
                        std::cout << ref << SEP << altOut << SEP << sampleCt << SEP << effectiveSampleCt << std::endl;
                    }
                    if (outfile) {
                        if (!forcingUnphased) {
                            writer->writeVariantSamples(*igdOutfile, pos, ref, alt, sampleList, isMissing, numCopies);
                            if (!vids.empty() && !noVariantIds) {
                                newVariantIds.emplace_back(vids.at(i));
                            }
                        } else {
                            std::vector<IGDSampleList> samplesByCopies(ploidy);
                            SampleT startSampleId = SAMPLE_INDEX_NOT_SET;
                            size_t numCopies = 0;
                            for (const auto& sampleId : sampleList) {
                                if (sampleId >= (startSampleId + ploidy) || startSampleId == SAMPLE_INDEX_NOT_SET) {
                                    if (numCopies > 0) {
                                        samplesByCopies[numCopies - 1].emplace_back(startSampleId / ploidy);
                                    }
                                    numCopies = 1;
                                    startSampleId = (sampleId / ploidy) * ploidy; // Round down by ploidy
                                } else {
                                    numCopies++;
                                }
                            }
                            if (numCopies > 0) {
                                samplesByCopies[numCopies - 1].emplace_back(startSampleId / ploidy);
                            }
                            for (uint8_t nc = 1; nc <= ploidy; nc++) {
                                const auto& copySampleList = samplesByCopies[nc - 1];
                                if (!copySampleList.empty()) {
                                    writer->writeVariantSamples(
                                        *igdOutfile, pos, ref, alt, copySampleList, isMissing, nc);
                                    if (!vids.empty() && !noVariantIds) {
                                        newVariantIds.emplace_back(vids.at(i));
                                    }
                                }
                            }
                        }
                    }
                    if (stats) {
                        variants++;
                        if (pos != lastPosition) {
                            sites++;
                            lastPosition = pos;
                        }
                        if (isMissing) {
                            missingRows++;
                        }
                        for (auto sampleId : sampleList) {
                            if (!isMissing) {
                                sampleToMuts.at(sampleId)++;
                            } else {
                                missingAlleles++;
                            }
                        }
                        samplesPerVariant.push_back(sampleCt);
                        sampleRefsTotal += sampleCt;
                        variantsByCopy.at(numCopies)++;
                    }
                }
            }

            if (stats) {
                std::cout << "Stats for " << *infile << std::endl;
                std::cout << "... in range " << bpStart << " - " << bpEnd << std::endl;
                std::cout << "  Variants in range: " << variants << std::endl;
                const double avgSamples = (double)sampleRefsTotal / (double)variants;
                std::cout << "  Average samples/var: " << avgSamples << std::endl;
                double stddevSamples = 0.0;
                for (auto count : samplesPerVariant) {
                    double diff = (double)count - avgSamples;
                    stddevSamples += (diff * diff);
                }
                stddevSamples = sqrt(stddevSamples / (double)variants);
                std::cout << "  Stddev samples/var: " << stddevSamples << std::endl;

                size_t mutRefsTotal = 0;
                for (size_t i = 0; i < sampleToMuts.size(); i++) {
                    mutRefsTotal += sampleToMuts[i];
                }
                const double avgMuts = (double)mutRefsTotal / (double)effectiveSampleCt;
                std::cout << "  Average var/sample: " << avgMuts << std::endl;
                double stddevMuts = 0.0;
                for (size_t i = 0; i < sampleToMuts.size(); i++) {
                    double diff = (double)sampleToMuts[i] - avgMuts;
                    stddevMuts += (diff * diff);
                }
                stddevMuts = sqrt(stddevMuts / (double)effectiveSampleCt);
                std::cout << "  Stddev var/sample: " << stddevMuts << std::endl;

                std::cout << "  Variants with missing data: " << missingRows << std::endl;
                std::cout << "  Total missing alleles: " << missingAlleles << std::endl;
                std::cout << "  Total unique sites: " << sites << std::endl;

                if (!igd.isPhased()) {
                    std::cout << "  Variants by num_copies:" << std::endl;
                    for (size_t i = 0; i < variantsByCopy.size(); i++) {
                        std::cout << "    num_copies=" << i << ", variants=" << variantsByCopy[i] << std::endl;
                    }
                }
            }
            if (outfile) {
                writer->writeIndex(*igdOutfile);
                writer->writeVariantInfo(*igdOutfile);
                std::vector<std::string> idsToUse;
                if (updateIndividualIds) {
                    idsToUse = readLines(*updateIndividualIds);
                } else {
                    const auto& indivIds = igd.getIndividualIds();
                    if (!indivIds.empty()) {
                        idsToUse.resize(numIndividuals);
                        for (size_t i = 0; i < indivIds.size(); i++) {
                            const size_t mappedIndex = mapIndividual.at(i);
                            if (mappedIndex != std::numeric_limits<size_t>::max()) {
                                idsToUse.at(mappedIndex) = indivIds[i];
                            }
                        }
                    }
                }
                if (!idsToUse.empty()) {
                    if (idsToUse.size() != numIndividuals) {
                        throw MalformedFile("Input file has invalid number of individual IDs");
                    }
                    writer->writeIndividualIds(*igdOutfile, idsToUse);
                }
                if (!newVariantIds.empty() && !noVariantIds) {
                    writer->writeVariantIds(*igdOutfile, newVariantIds);
                }
                igdOutfile->seekp(0);
                writer->writeHeader(*igdOutfile, *infile, igd.getDescription());
            }
        }
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}
