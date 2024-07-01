/* Convert VCF to IIGT.
 *
 * Usage:
 *  vcfconv <vcf-file> <output-file> -copy-ids
 */
#include <iostream>

#include "picovcf.hpp"

using namespace picovcf;

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: vcfconv <vcf-file> <output-file> [-copy-ids]" << std::endl;
        return 1;
    }

    bool emitIndividualIds = false;
    bool emitVariantIds = false;
    const std::string infile(argv[1]);
    const std::string outfile(argv[2]);
    if (argc > 3) {
        std::string arg3 = argv[3];
        if (arg3 == "-copy-ids") {
            emitIndividualIds = true;
            emitVariantIds = true;
        } else {
            std::cerr << "Unrecognized flag \"" << arg3 << "\"" << std::endl;
            return 1;
        }
    }
    vcfToIGD(infile, outfile, "", true, emitIndividualIds, emitVariantIds);
    return 0;
}