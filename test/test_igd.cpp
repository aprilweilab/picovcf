#include <gtest/gtest.h>

#include "picovcf.hpp"

#include <string>

using namespace picovcf;

extern const std::string getMISSING_DATA_EXAMPLE_FILE();

TEST(IGD, MissingData) {
    std::string igdFileName = "missing_data.igd";
    vcfToIGD(getMISSING_DATA_EXAMPLE_FILE(), igdFileName);

    IGDData igdFile(igdFileName);

    ASSERT_EQ(igdFile.numVariants(), 12); // 2 of the variants have multiple alt alleles
    size_t numMissing = 0;
    for (size_t i = 0; i < igdFile.numVariants(); i++) {
        bool isMissing = false;
        igdFile.getPosition(i, isMissing);
        auto sampleList = igdFile.getSamplesWithAlt(i);
        if (isMissing) {
            numMissing++;
        }
    }
    // Every variant is missing
    ASSERT_EQ(numMissing, 5);
}
