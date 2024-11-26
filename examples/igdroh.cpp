/* Example: Compute ROH from an unphased IGD file.
 *
 * Usage:
 *  igdroh <unphased IGD file>
 */
#include <iostream>

#include "picovcf.hpp"

using namespace picovcf;

constexpr size_t THRESHOLD = 500000; // 500kbp

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Please pass in an input file" << std::endl;
        return 1;
    }

    const std::string filename(argv[1]);
    IGDData igd(filename);
    if (igd.isPhased()) {
        std::cerr << "This example only works with unphased inputs" << std::endl;
        return 2;
    }

    // Just emits ROH as spans (individual, start_bp, end_bp) as an example. Only emits
    // runs that are longer than THRESHOLD above.
    std::vector<uint64_t> lastHeteroSitePerIndiv(igd.numIndividuals(), 0);
    static constexpr char SEP = '\t';
    std::cout << "INDIV" << SEP << "START_BP" << SEP << "END_BP" << std::endl;
    uint64_t position = 0;
    // We iterate all variants, but skip the homozygotes. The span between two heterozygotes
    // is the length of the ROH, and if it exceeds our threshold we emit it.
    for (size_t i = 0; i < igd.numVariants(); i++) {
        bool isMissing = false;
        uint8_t numCopies = 0;
        position = igd.getPosition(i, isMissing, numCopies);
        assert(numCopies > 0); // This should always be true for unphased data
        if (numCopies == 1) {
            for (const auto indivIndex : igd.getSamplesWithAlt(i)) {
                const uint64_t homozygSpan = (position - lastHeteroSitePerIndiv[indivIndex]);
                if (homozygSpan >= THRESHOLD) {
                    std::cout << indivIndex << SEP << lastHeteroSitePerIndiv[indivIndex] + 1 << SEP << position - 1
                              << std::endl;
                }
                lastHeteroSitePerIndiv[indivIndex] = position;
            }
        }
    }
    // Emit any ROH segments that go to the end of the genome range.
    for (size_t indivIndex = 0; indivIndex < igd.numIndividuals(); indivIndex++) {
        const uint64_t homozygSpan = (position - lastHeteroSitePerIndiv[indivIndex]);
        if (homozygSpan >= THRESHOLD) {
            std::cout << indivIndex << SEP << lastHeteroSitePerIndiv[indivIndex] + 1 << SEP << position - 1
                      << std::endl;
        }
        lastHeteroSitePerIndiv[indivIndex] = position;
    }

    return 0;
}