#include <gtest/gtest.h>

#include "picovcf.hpp"

#include <string>

using namespace picovcf;

extern const std::string getMISSING_DATA_EXAMPLE_FILE();

TEST(IGD, MissingData) {
    std::string igdFileName = "missing_data.igd";
    vcfToIGD(getMISSING_DATA_EXAMPLE_FILE(), igdFileName);

    IGDData igdFile(igdFileName);

    ASSERT_EQ(igdFile.numVariants(), 7); // 2 of the variants have multiple alt alleles
    const auto& missingData = igdFile.getMissingData();
    // Every variant is missing
    ASSERT_EQ(missingData.size(), 5);
}
