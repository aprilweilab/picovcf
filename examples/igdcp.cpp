/* Copy from one IGD file to another.
 *
 * Usage:
 *  igdcp <infile> <outfile>
 */
#include <iostream>
#include <cmath>
#include <iomanip>

#include "picovcf.hpp"

using namespace picovcf;

int main(int argc, char *argv[]) {
    std::cout << std::fixed << std::setprecision(4);
    if (argc < 3) {
        std::cerr << "Usage: igdcp <infile> <outfile>" << std::endl;
        return 1;
    }

    const std::string infile(argv[1]);
    const std::string outfile(argv[2]);

    IGDData igd(infile);
    uint64_t numIndividuals = igd.numIndividuals();
    uint64_t ploidy = igd.getPloidy();

    // Example of copying data from one IGD to another. Copy all the variants and then write the
    // index and optional metadata at the end. We write the header twice because it contains file
    // pointers that get updated by the other "writeXXX()" methods.
    std::ofstream igdOutfile(outfile, std::ios::binary);
    IGDWriter writer(ploidy, igd.numIndividuals(), igd.isPhased());
    writer.writeHeader(igdOutfile, infile, igd.getDescription());
    for (size_t i = 0; i < igd.numVariants(); i++) {
        bool isMissing = false;
        const auto position = igd.getPosition(i, isMissing);
        writer.writeVariantSamples(igdOutfile,
                                   position,
                                   igd.getRefAllele(i),
                                   igd.getAltAllele(i),
                                   igd.getSamplesWithAlt(i),
                                   isMissing);
    }
    writer.writeIndex(igdOutfile);
    writer.writeVariantInfo(igdOutfile);
    writer.writeIndividualIds(igdOutfile, igd.getIndividualIds());
    writer.writeVariantIds(igdOutfile, igd.getVariantIds());
    igdOutfile.seekp(0);
    writer.writeHeader(igdOutfile, infile, igd.getDescription());

    return 0;
}
