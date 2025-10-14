/* Given an unphased IGD, "phase" it by assigning all heterozygotes to haplotype copy 1,
 * except for the ones listed in another IGD file (which are put on hap copy 2)
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
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <sys/stat.h>

#include "picovcf.hpp"

using namespace picovcf;

const std::string DIR_SEP = "/";

struct SampleData {
    IGDSampleList hetSamples;
    IGDSampleList homSamples;
    IGDSampleList missing;
    std::string ref;
    std::string alt;
    size_t position;
};

static bool cmpSampleData(const SampleData& sd1, const SampleData& sd2) {
    if (sd1.ref == sd2.ref) {
        return sd1.alt < sd2.alt;
    }
    return sd1.ref < sd2.ref;
}

inline void setOrAppendSamples(IGDSampleList& target, IGDSampleList source) {
    if (target.empty()) {
        target = std::move(source);
    } else {
        for (auto sampleIdx : source) {
            target.push_back(sampleIdx);
        }
        std::sort(target.begin(), target.end());
    }
}

// Return all sample data for the current site, grouped by alleles. The resulting vector is in ascending order
// of (ref, alt).
std::vector<SampleData> collectSiteByAlleles(IGDData& igd, size_t& variantIndex) {
    PICOVCF_RELEASE_ASSERT(variantIndex < igd.numVariants());
    std::vector<SampleData> result;
    std::vector<size_t> siteVariants = igd.collectNextSite(variantIndex);
    for (size_t s : siteVariants) {
        uint8_t numCopies = 0;
        bool isMissing = false;
        auto pos = igd.getPosition(s, isMissing, numCopies);
        const auto& ref = igd.getRefAllele(s);
        const auto& alt = igd.getAltAllele(s);
        size_t resultIdx = 0;
        for (resultIdx = 0; resultIdx < result.size(); resultIdx++) {
            if (result[resultIdx].ref == ref && result[resultIdx].alt == alt) {
                break;
            }
        }
        if (resultIdx == result.size()) {
            result.push_back({{}, {}, {}, ref, alt, pos});
        }
        if (numCopies == 2) {
            setOrAppendSamples(result[resultIdx].homSamples, std::move(igd.getSamplesWithAlt(s)));
        } else {
            PICOVCF_RELEASE_ASSERT(numCopies < 2);
            if (isMissing) {
                setOrAppendSamples(result[resultIdx].missing, std::move(igd.getSamplesWithAlt(s)));
            } else {
                setOrAppendSamples(result[resultIdx].hetSamples, std::move(igd.getSamplesWithAlt(s)));
            }
        }
    }
    variantIndex = siteVariants.back() + 1;
    std::sort(result.begin(), result.end(), cmpSampleData);
    return std::move(result);
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: igd_flip <unphased IGD> <flip IGD> <output IGD>" << std::endl;
        return 1;
    }

    // Not a VCF, then assume it is IGD and load the header.
    std::string infile = argv[1];
    IGDData unphasedIgd(infile);
    IGDData flipIgd(argv[2]);
    std::ofstream outStream(argv[3], std::ios::binary);
    IGDWriter writer(2, unphasedIgd.numIndividuals(), true);
    writer.writeHeader(outStream, infile, unphasedIgd.getDescription());

    size_t flipCount = 0;
    size_t unphIdx = 0;
    size_t flipIdx = 0;
    while (unphIdx < unphasedIgd.numVariants()) {
        std::vector<SampleData> data = collectSiteByAlleles(unphasedIgd, unphIdx);
        PICOVCF_RELEASE_ASSERT(!data.empty());
        const size_t position = data.front().position;

        bool flipMissing = false;
        size_t flipPos = INVALID_POSITION;
        // Catch up to our current position if needed.
        while (flipIdx < flipIgd.numVariants() && (flipPos = flipIgd.getPosition(flipIdx, flipMissing)) < position) {
            PICOVCF_RELEASE_ASSERT(!flipMissing);
            std::cout << "Skipping flip idx " << flipIdx << "\n";
            flipIdx++;
        }

        std::vector<SampleData> flipData =
            (flipIdx < flipIgd.numVariants() && flipPos == position) ? collectSiteByAlleles(flipIgd, flipIdx) : std::vector<SampleData>();
        for (auto& sampleData : data) {
            // Handle any missing data first.
            if (!sampleData.missing.empty()) {
                writer.writeVariantSamples(
                    outStream, sampleData.position, sampleData.ref, sampleData.alt, sampleData.missing, true);
            }

            // Convert the homozygotes into haplotypes.
            IGDSampleList haplotypeSamples;
            for (const auto& indivIdx : sampleData.homSamples) {
                haplotypeSamples.emplace_back(2 * indivIdx);
                haplotypeSamples.emplace_back((2 * indivIdx) + 1);
            }

            // Now look to see if there is any flip data.
            size_t j = 0;
            for (j = 0; j < flipData.size(); j++) {
                if (flipData[j].ref == sampleData.ref && flipData[j].alt == sampleData.alt) {
                    break;
                }
            }
            // There is. Do the flipping in place against sampleData.
            size_t x = 0;
            if (j < flipData.size()) {
                flipCount++;
                PICOVCF_RELEASE_ASSERT(flipData[j].homSamples.empty());
                PICOVCF_RELEASE_ASSERT(flipData[j].missing.empty());
                const auto& toFlip = flipData[j].hetSamples;
                size_t y = 0;
                while (x < sampleData.hetSamples.size() && y < toFlip.size()) {
                    if (toFlip[y] == sampleData.hetSamples[x]) {
                        haplotypeSamples.emplace_back((sampleData.hetSamples[x] * 2) + 1); // Right hap
                        x++;
                        y++;
                    } else if (toFlip[y] < sampleData.hetSamples[x]) {
                        y++;
                    } else {
                        haplotypeSamples.emplace_back((sampleData.hetSamples[x] * 2) + 0); // Left hap
                        x++;
                    }
                }
            }
            // Cleanup remainder
            while (x < sampleData.hetSamples.size()) {
                haplotypeSamples.emplace_back((sampleData.hetSamples[x] * 2) + 0);
                x++;
            }
            std::sort(haplotypeSamples.begin(), haplotypeSamples.end());
            writer.writeVariantSamples(
                outStream, sampleData.position, sampleData.ref, sampleData.alt, haplotypeSamples, false);
        }
    }
    std::cout << "Flipped: " << flipCount << "\n";
    std::cout << "Should have flipped: " << flipIgd.numVariants() << "\n";

    writer.writeIndex(outStream);
    writer.writeVariantInfo(outStream);
    writer.writeIndividualIds(outStream, unphasedIgd.getIndividualIds());
    // TODO: write variant IDs
    // writer->writeVariantIds(*igdOutfile, newVariantIds);
    outStream.seekp(0);
    writer.writeHeader(outStream, infile, unphasedIgd.getDescription());

    return 0;
}
