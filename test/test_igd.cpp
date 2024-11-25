#include <gtest/gtest.h>

#include "picovcf.hpp"

#include <string>

using namespace picovcf;

extern const std::string getMISSING_DATA_EXAMPLE_FILE();
extern const std::string getUNPHASED_DATA_EXAMPLE_FILE();

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

TEST(IGD, UnphasedData) {
    std::string igdFileName = "unphased_data.igd";
    vcfToIGD(getUNPHASED_DATA_EXAMPLE_FILE(), igdFileName);

    IGDData igdFile(igdFileName);

    ASSERT_EQ(igdFile.numVariants(), 7);
    size_t totalOnes = 0;
    size_t totalTwos = 0;
    for (size_t i = 0; i < igdFile.numVariants(); i++) {
        bool isMissing = false;
        uint8_t numCopies = 0;
        igdFile.getPosition(i, isMissing, numCopies);
        ASSERT_GT(numCopies, 0);
        ASSERT_LE(numCopies, 2);
        ASSERT_FALSE(isMissing);
        auto sampleList = igdFile.getSamplesWithAlt(i);
        for (const auto index : sampleList) {
            ASSERT_LT(index, igdFile.numIndividuals());
        }
        if (numCopies == 1) {
            totalOnes += sampleList.size();
        } else if (numCopies == 2) {
            totalTwos += sampleList.size();
        }
    }
    ASSERT_EQ(totalOnes, 90+325+146+129);
    ASSERT_EQ(totalTwos, 2+2+1);
}
