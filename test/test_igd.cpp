#include <gtest/gtest.h>

#include "picovcf.hpp"

#include <string>

using namespace picovcf;

extern const std::string getMISSING_DATA_EXAMPLE_FILE();
extern const std::string getUNPHASED_DATA_EXAMPLE_FILE();
extern const std::string getMIXED_DATA_EXAMPLE_FILE();
extern const std::string getMSPRIME_EXAMPLE_FILE();

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

    // 7 + the 3 missingness variants
    ASSERT_EQ(igdFile.numVariants(), 10);
    size_t totalOnes = 0;
    size_t totalTwos = 0;
    size_t totalMiss = 0;
    for (size_t i = 0; i < igdFile.numVariants(); i++) {
        bool isMissing = false;
        uint8_t numCopies = 0;
        igdFile.getPosition(i, isMissing, numCopies);
        ASSERT_TRUE(isMissing || numCopies > 0);
        ASSERT_LE(numCopies, 2);
        auto sampleList = igdFile.getSamplesWithAlt(i);
        for (const auto index : sampleList) {
            ASSERT_LT(index, igdFile.numIndividuals());
        }
        if (isMissing) {
            totalMiss += sampleList.size();
        } else if (numCopies == 1) {
            totalOnes += sampleList.size();
        } else if (numCopies == 2) {
            totalTwos += sampleList.size();
        }
    }
    ASSERT_EQ(totalOnes, 90+325+146+129);
    ASSERT_EQ(totalTwos, 2+2+1);
    ASSERT_EQ(totalMiss, 4);
}

TEST(IGD, MixedPloidyData) {
    std::string igdFileName = "mixed_ploidy.igd";
    vcfToIGD(getMIXED_DATA_EXAMPLE_FILE(), igdFileName, "", false, true, true, false, 2);

    IGDData igdFile(igdFileName);

	const size_t numSamples = igdFile.numSamples();
    ASSERT_EQ(numSamples, 20000);
    ASSERT_EQ(igdFile.numVariants(), 4);
    // First 2 variants were haploid data in the VCF file, so they should have been "mirrored"
    // to diploid (so both alleles are the same).
    for (size_t i = 0; i < 2; i++) {
        bool isMissing = false;
        uint8_t numCopies = 0;
        igdFile.getPosition(i, isMissing, numCopies);
        ASSERT_FALSE(isMissing);
        ASSERT_EQ(numCopies, 0);
        auto sampleList = igdFile.getSamplesWithAlt(i);
        SampleT prevSample = SAMPLE_INDEX_NOT_SET;
        for (const auto index : sampleList) {
            if (index % 2 == 0) {
                prevSample = index;
            } else {
                ASSERT_EQ(prevSample, index - 1);
            }
        }
    }
}

TEST(IGD, LowerBound) {
    std::string igdFileName = "lower_bound.igd";
    vcfToIGD(getMSPRIME_EXAMPLE_FILE(), igdFileName);

    IGDData igdFile(igdFileName);

    ASSERT_EQ(igdFile.numVariants(), 7);
    size_t variantIndex = igdFile.lowerBoundPosition(0);
    ASSERT_EQ(variantIndex, 0);
    variantIndex = igdFile.lowerBoundPosition(56554);
    ASSERT_EQ(variantIndex, 3);
    ASSERT_EQ(igdFile.getPosition(variantIndex), 56554);
    variantIndex = igdFile.lowerBoundPosition(56553);
    ASSERT_EQ(variantIndex, 3);
    variantIndex = igdFile.lowerBoundPosition(56812);
    ASSERT_EQ(variantIndex, 5);
    variantIndex = igdFile.lowerBoundPosition(100000000);
    ASSERT_EQ(variantIndex, igdFile.numVariants());
}

