/* Tools for converting to and manipulating IGD files.
 *
 * Run "igdtools --help" for usage.
 */
#include <iostream>
#include <cmath>
#include <iomanip>
#include <limits>

#include "third-party/args.hxx"
#include "picovcf.hpp"

using namespace picovcf;

inline void emitAllele(VariantT alleleIndex, std::ostream& out) {
    if (alleleIndex == MISSING_VALUE) {
        out << "? ";
    } else {
        out << alleleIndex << " ";
    }
}

inline bool ends_with(std::string const &string1, std::string const &string2) {
    if (string1.length() < string2.length()) {
        return false;
    }
    return (string2 == string1.substr(string1.length() - string2.length()));
}

template <typename Out>
inline void split(const std::string &s, char delim, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

inline std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return std::move(elems);
}

int main(int argc, char *argv[]) {
    std::cout << std::fixed << std::setprecision(4);

    args::ArgumentParser parser("Process or create IGD files.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::Positional<std::string> infile(parser, "input_file", "The input file (.vcf, .vcf.gz, or .igd)");
    args::ValueFlag<std::string> outfile(
        parser, "output", "The output file to produce.", {'o', "out"});
    args::ValueFlag<std::string> outputDescription(
        parser, "outputDescription", "The description string to include in the IGD output.", {"description"});
    args::ValueFlag<std::string> range(
        parser, "range", "Restrict to the given base-pair range (inclusive).", {'r', "range"});
    args::ValueFlag<std::string> frange(
        parser, "frange", "Restrict to variants with frequency in the range (inclusive, exclusive).", {'f', "frange"});
    args::Flag info(
        parser, "info", "Display information from the IGD header.", {'i', "info"});
    args::Flag individuals(
        parser, "individuals", "Emit the mapping from individual index to ID.", {"individuals"});
    args::Flag variants(
        parser, "variants", "Emit the mapping from variant index to ID.", {"variants"});
    args::Flag stats(
        parser, "stats", "Emit some simple statistics about the distribution of samples, sites, and variants.", {'s', "stats"});
    args::Flag alleles(
        parser, "alleles", "Emit allele frequencies.", {'a', "alleles"});
    args::Flag noIndividualIds(
        parser, "noIndividualIds", "Do not emit IDs for individuals in the resulting IGD file.", {"no-indiv-ids"}
    );
    args::Flag noVariantIds(
        parser, "noVariantIds", "Do not emit IDs for variants in the resulting IGD file.", {"no-var-ids"}
    );
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

#define UNSUPPORTED_FOR_VCF(parameter, parameterName) do { \
    if ((bool)(parameter)) { \
        std::cerr << "Parameter " << (parameterName) << " is not supported for VCF file conversion" << std::endl; \
        return 1; \
    } \
} while(0)

    std::string description = outputDescription ? *outputDescription : "";

    const bool isVcf = ends_with(*infile, ".vcf") || ends_with(*infile, ".vcf.gz");
    if (isVcf) {
        if (!outfile) {
            std::cerr << "VCF input is only supported for conversion to IGD. Use --output." << std::endl;
            return 1;
        }
        UNSUPPORTED_FOR_VCF(range, "--range");

        const bool emitIndividualIds = !noIndividualIds;
        const bool emitVariantIds = !noVariantIds;
        vcfToIGD(*infile, *outfile, description, true, emitIndividualIds, emitVariantIds);
        return 0;
    }

    // Not a VCF, then assume it is IGD and load the header.
    IGDData igd(*infile);

    size_t bpStart = 0;
    size_t bpEnd = std::numeric_limits<size_t>::max();
    if (range) {
        const char *rangec = range->c_str();
        const char *rangecEnd = range->c_str() + range->size();
        char *endptr = nullptr;
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

    double fLower = 0.0;
    double fUpper = 1.1;
    if (frange) {
        const char *frangec = frange->c_str();
        const char *frangecEnd = frange->c_str() + frange->size();
        char *endptr = nullptr;
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

    const size_t numSamples = igd.numSamples();

    if (info) {
        std::cout << "Header information for " << *infile << std::endl;
        std::cout << "  Variants: " << igd.numVariants() << std::endl;
        std::cout << "  Individuals: " << igd.numIndividuals() << std::endl;
        std::cout << "  Ploidy: " << igd.getPloidy() << std::endl;
        std::cout << "  Phased?: " << (igd.isPhased() ? "true" : "false") << std::endl;
        std::cout << "  Source: " << igd.getSource() << std::endl;
        std::cout << "  Genome range: " << igd.getPosition(0)
                                        << "-" << igd.getPosition(igd.numVariants()-1) << std::endl;
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

#define CONDITION_PRINT(condition, message) do { \
    if (condition) { \
        std::cout << message; \
    } \
} while(0)

    std::shared_ptr<std::ofstream> igdOutfile;
    std::shared_ptr<IGDWriter> writer;
    if (outfile) {
        uint64_t numIndividuals = igd.numIndividuals();
        uint64_t ploidy = igd.getPloidy();

        igdOutfile = std::make_shared<std::ofstream>(*outfile, std::ios::binary);
        writer = std::make_shared<IGDWriter>(ploidy, igd.numIndividuals(), igd.isPhased());
        writer->writeHeader(*igdOutfile, *infile, igd.getDescription());
    }

    const bool iterateSamples = (bool)stats || (bool)alleles || (bool)outfile;
    if (iterateSamples) {
        static constexpr char SEP = '\t';
        CONDITION_PRINT(alleles, "POSITION" << SEP << "REF" << SEP << "ALT" << SEP << "ALT COUNT" << SEP << "TOTAL" << std::endl);

        size_t sites = 0;
        bool _ignore = false;
        size_t variants = 0;
        size_t missingRows = 0;
        size_t missingAlleles = 0;
        std::vector<size_t> samplesPerVariant;
        size_t sampleRefsTotal = 0;
        std::vector<size_t> sampleToMuts(igd.numSamples()); // Counts muts per sample
        
        size_t lastPosition = std::numeric_limits<size_t>::max();
        for (size_t i = 0; i < igd.numVariants(); i++) {
            bool isMissing = false;
            auto pos = igd.getPosition(i, isMissing);
            if (pos >= bpStart && pos <= bpEnd) {
                auto sampleList = igd.getSamplesWithAlt(i);
                const auto sampleCt = sampleList.size();
                const double freq = (double)sampleCt/(double)numSamples;
                if (freq < fLower || freq >= fUpper) {
                    continue;
                }
                const auto& ref = igd.getRefAllele(i);
                const auto& alt = igd.getAltAllele(i);
                CONDITION_PRINT(alleles, pos << SEP << ref << SEP << alt << SEP << sampleCt << SEP << igd.numSamples() << std::endl);
                if (outfile) {
                    writer->writeVariantSamples(*igdOutfile,
                                                pos,
                                                ref,
                                                alt,
                                                sampleList,
                                                isMissing);
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
            stddevSamples = sqrt(stddevSamples/(double)variants);
            std::cout << "  Stddev samples/var: " << stddevSamples << std::endl;

            size_t mutRefsTotal = 0;
            for (size_t i = 0; i < sampleToMuts.size(); i++) {
                mutRefsTotal += sampleToMuts[i];
            }
            const double avgMuts = (double)mutRefsTotal / (double)igd.numSamples();
            std::cout << "  Average var/sample: " << avgMuts << std::endl;
            double stddevMuts = 0.0;
            for (size_t i = 0; i < sampleToMuts.size(); i++) {
                double diff = (double)sampleToMuts[i] - avgMuts;
                stddevMuts += (diff*diff);
            }
            stddevMuts = sqrt(stddevMuts/(double)igd.numSamples());
            std::cout << "  Stddev var/sample: " << stddevMuts << std::endl;
        
            std::cout << "  Variants with missing data: " << missingRows << std::endl;
            std::cout << "  Total missing alleles: " << missingAlleles << std::endl;
            std::cout << "  Total unique sites: " << sites << std::endl;
        }
        if (outfile) {
            writer->writeIndex(*igdOutfile);
            writer->writeVariantInfo(*igdOutfile);
            writer->writeIndividualIds(*igdOutfile, igd.getIndividualIds());
            writer->writeVariantIds(*igdOutfile, igd.getVariantIds());
            igdOutfile->seekp(0);
            writer->writeHeader(*igdOutfile, *infile, igd.getDescription());
        }
    }

    return 0;
}
