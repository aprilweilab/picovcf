/* VCF[.gz] pretty-printing tool. This mostly just exists for demonstrating and
 * testing the usage of picovcf.
 *
 * Usage:
 *  vcfpp <command> <file>
 *
 * where command is one of ["stats", "matrix"]
 */
#include <iostream>

#include "picovcf.hpp"

using namespace picovcf;

inline void emitAllele(VariantT alleleIndex, std::ostream& out) {
    if (alleleIndex == MISSING_VALUE) {
        out << "? ";
    } else {
        out << alleleIndex << " ";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Please pass in a command and an input file (in that order)" << std::endl;
        return 1;
    }

    const std::string command(argv[1]);
    const std::string filename(argv[2]);

    VCFFile vcf(filename);

    if (command == "fast_stats") {
        // Count numVariants() is the slow part, so we leave that off!
        std::cout << "Stats for " << filename << std::endl;
        std::cout << "  Has index? " << (vcf.isUsingIndex() ? "yes" : "no") << std::endl;
        std::cout << "  Individuals: " << vcf.numIndividuals() << std::endl;
        std::cout << "  Version: " << vcf.getMetaInfo(VCFFile::META_FILE_FORMAT) << std::endl;
        std::cout << "  Source: " << vcf.getMetaInfo("source") << std::endl;
        std::cout << "  Genome range: " << vcf.getGenomeRange().first << "-" << vcf.getGenomeRange().second
                  << std::endl;
        if (vcf.nextVariant()) {
            std::cout << "  Has genotype data? " << (vcf.currentVariant().hasGenotypeData() ? "yes" : "no")
                      << std::endl;
        }
    } else if (command == "stats") {
        std::cout << "Stats for " << filename << std::endl;
        std::cout << "  Has index? " << (vcf.isUsingIndex() ? "yes" : "no") << std::endl;
        std::cout << "  Variants: " << vcf.numVariants() << std::endl;
        std::cout << "  Individuals: " << vcf.numIndividuals() << std::endl;
        std::cout << "  Version: " << vcf.getMetaInfo(VCFFile::META_FILE_FORMAT) << std::endl;
        std::cout << "  Source: " << vcf.getMetaInfo("source") << std::endl;
        std::cout << "  Genome range: " << vcf.getGenomeRange().first << "-" << vcf.getGenomeRange().second
                  << std::endl;
        if (vcf.nextVariant()) {
            std::cout << "  Has genotype data? " << (vcf.currentVariant().hasGenotypeData() ? "yes" : "no")
                      << std::endl;
        }
    } else if (command == "matrix") {
        std::cout << "Rows: variants" << std::endl;
        std::cout << "Columns: samples" << std::endl;
        // This assumes/checks for phased data.
        vcf.seekBeforeVariants();
        while (vcf.nextVariant()) {
            VCFVariantView& variant = vcf.currentVariant();
            std::vector<AlleleT> gt = variant.getGenotypeArray();
            for (const auto allele : gt) {
                emitAllele(allele, std::cout);
            }
            std::cout << std::endl;
        }
    }
    return 0;
}
