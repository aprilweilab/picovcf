/* Convert VCF to IIGT.
 *
 * Usage:
 *  vcfconv <vcf-file> <output-file>
 */
#include <iostream>

#include "picovcf.hpp"

using namespace picovcf;

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Please pass in an input file and output file name" << std::endl;
        return 1;
    }

    const std::string infile(argv[1]);
    const std::string outfile(argv[2]);
    vcfToIGD(infile, outfile);
    return 0;
}