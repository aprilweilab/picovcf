/* IGD pretty-printing tool. This mostly just exists for demonstrating and testing
 * the usage of picovcf.
 *
 * Usage:
 *  igdpp <command> <file>
 *
 * where command is one of ["stats", "range_start"]
 */
#include <iostream>
#include <cmath>
#include <iomanip>

#include "picovcf.hpp"

using namespace picovcf;

inline void emitAllele(VariantT alleleIndex, std::ostream& out) {
    if (alleleIndex == MISSING_VALUE) {
        out << "? ";
    } else {
        out << alleleIndex << " ";
    }
}

int main(int argc, char *argv[]) {
    std::cout << std::fixed << std::setprecision(4);
    if (argc < 3) {
        std::cerr << "Please pass in a command and an input file (in that order)" << std::endl;
        return 1;
    }

    const std::string command(argv[1]);
    const std::string filename(argv[2]);

    IGDData igd(filename);

    if (command == "stats") {
        std::cout << "Stats for " << filename << std::endl;
        std::cout << "  Variants: " << igd.numVariants() << std::endl;
        std::cout << "  Individuals: " << igd.numIndividuals() << std::endl;
        std::cout << "  Ploidy: " << igd.getPloidy() << std::endl;
        std::cout << "  Phased?: " << (igd.isPhased() ? "true" : "false") << std::endl;
        std::cout << "  Source: " << igd.getSource() << std::endl;
        std::cout << "  Genome range: " << igd.getPosition(0)
                                        << "-" << igd.getPosition(igd.numVariants()-1) << std::endl;
        std::cout << "  Has individual IDs? " << (igd.getIndividualIds().empty() ? "No" : "Yes") << std::endl;
    } else if (command == "range_stats") {
        std::cout << "Stats for " << filename << std::endl;
        bool _ignore = false;
        size_t firstPos = igd.getPosition(0, _ignore);
        size_t lastPos = igd.getPosition(igd.numVariants()-1, _ignore);
        uint32_t start = (argc > 3) ? atoi(argv[3]) : firstPos;
        uint32_t end = (argc > 4) ? atoi(argv[4]) : lastPos + 1;
        size_t variants = 0;
        size_t missingRows = 0;
        size_t missingAlleles = 0;
        std::vector<size_t> samplesPerVariant;
        size_t sampleRefsTotal = 0;
        std::vector<size_t> sampleToMuts(igd.numSamples()); // Counts muts per sample
        for (size_t i = 0; i < igd.numVariants(); i++) {
            bool isMissing = false;
            auto pos = igd.getPosition(i, isMissing);
            if (pos >= start && pos < end) {
                variants++;
                auto sampleList = igd.getSamplesWithAlt(i);
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
                auto sampleCt = sampleList.size();
                samplesPerVariant.push_back(sampleCt);
                sampleRefsTotal += sampleCt;
            }
        }
        std::cout << "... in range " << start << " - " << end << std::endl;
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


    } else if (command == "freq") {
        static constexpr char SEP = '\t';
        std::cout << "POSITION" << SEP << "REF" << SEP << "ALT" << SEP << "ALT COUNT" << SEP << "TOTAL" << std::endl;
        for (size_t i = 0; i < igd.numVariants(); i++) {
            bool isMissing = false;
            auto pos = igd.getPosition(i, isMissing);
            auto sampleList = igd.getSamplesWithAlt(i);
            std::cout << pos << SEP << igd.getRefAllele(i) << SEP
                << igd.getAltAllele(i) << SEP << sampleList.size() << SEP << igd.numSamples() << std::endl;
        }
    }
    return 0;
}
