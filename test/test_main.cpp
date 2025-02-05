#include <gtest/gtest.h>

std::string getExampleDir() {
    const char * exampleDir = getenv("EXAMPLE_VCFS");
    if (nullptr == exampleDir) {
        return std::string(".");
    }
    return std::string(exampleDir);
}

const std::string getVCF42_EXAMPLE_FILE() {
    return getExampleDir() + "/vcf42.example.vcf";
}

const std::string getMSPRIME_EXAMPLE_FILE() {
    return getExampleDir() + "/msprime.example.vcf";
}

const std::string getMSPRIMEGZ_EXAMPLE_FILE() {
    return getExampleDir() + "/msprime.example.vcf.gz";
}

const std::string getMISSING_DATA_EXAMPLE_FILE() {
    return getExampleDir() + "/missing.vcf";
}

const std::string getUNPHASED_DATA_EXAMPLE_FILE() {
    return getExampleDir() + "/unphased.example.vcf";
}

const std::string getHAPLOID_DATA_EXAMPLE_FILE() {
    return getExampleDir() + "/haploid.example.vcf";
}

const std::string getMIXED_DATA_EXAMPLE_FILE() {
    return getExampleDir() + "/mixploidy.example.vcf";
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
