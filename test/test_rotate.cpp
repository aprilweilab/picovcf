#include <gtest/gtest.h>

#include "picovcf.hpp"

#include <string>

extern const std::string getVCF42_EXAMPLE_FILE();

using namespace picovcf;
using Matrix = std::vector<std::vector<size_t>>;
    
Matrix getMatrix(const VCFRotatedWindow &rotation) {
    Matrix rv(rotation.sampleToMutation.size());
    for (size_t i = 0; i < rotation.sampleToMutation.size(); i++) {
        rv[i].resize(rotation.parsedVariants.size());
        for (auto mutation : rotation.sampleToMutation[i]) {
            rv[i][mutation.first] = mutation.second;
        }
    }
    return rv;
}

Matrix combineColumns(const Matrix& a, const Matrix& b) {
    Matrix result(b.size());
    for (size_t i = 0; i < b.size(); i++) {
        if (!a.empty()) {
            for (auto v : a[i]) {
                result[i].push_back(v);
            }
        }
        for (auto v : b[i]) {
            result[i].push_back(v);
        }
    }
    return result;
}

Matrix combineRows(const Matrix& a, const Matrix& b) {
    Matrix result(a.size() + b.size());
    for (size_t i = 0; i < a.size(); i++) {
        for (auto v : a[i]) {
            result[i].push_back(v);
        }
    }
    for (size_t i = 0; i < b.size(); i++) {
        for (auto v : b[i]) {
            result[a.size()+i].push_back(v);
        }
    }
    return result;
}

TEST(Rotate, SpecExample) {
    // This example is 5 (rows) x 3 (columns). We want to rotate it to be 3 x 5, conceptually.
    VCFFile vcf(getVCF42_EXAMPLE_FILE());
    
    const Matrix expectedRotation = Matrix({
        {0, 0, 1, 0, 0},
        {0, 0, 2, 0, 1},
        {1, 0, 2, 0, 0},
        {0, 1, 1, 0, 2},
        {1, 0, 2, 0, 1},
        {1, 0, 2, 0, 1},
    });
    
    // Do it all at once: just rotate the whole file since it is small.
    VCFRotatedWindow rotation;
    vcf.getRotatedWindow({0, 3},       // Individual range
                         {0, 1234568}, // Genome range
                         rotation);
    ASSERT_TRUE(rotation.isDiploid);
    ASSERT_EQ(rotation.firstIndividual, 0);
    std::vector<size_t> positions;
    for (auto v : rotation.parsedVariants) {
        positions.push_back(v.position);
    }
    ASSERT_EQ(positions, std::vector<size_t>({14370, 17330, 1110696, 1230237, 1234567}));
    Matrix rotatedMatrix = getMatrix(rotation);
    ASSERT_EQ(rotatedMatrix, expectedRotation);

    // Now do it one column at a time
    Matrix windowedCols;
    for (size_t i = 0; i < 3; i++) {
        VCFRotatedWindow rotation;
        vcf.getRotatedWindow({i, i+1},      // Individual range
                             {0, 1234568}, // Genome range
                             rotation);
        Matrix partialMatrix = getMatrix(rotation);
        windowedCols = combineRows(windowedCols, partialMatrix);
    }
    ASSERT_EQ(windowedCols, expectedRotation);

    // Now do it one row at a time
    Matrix windowedRows;
    std::vector<RangePair> genomeRanges = {
        {0, 14371},
        {14371, 1110696},
        {1110696, 1200000},
        {1200000, 1234567},
        {1234567, 1234568},
    };
    FileOffset lastPosition = {0, 0};
    for (size_t i = 0; i < 5; i++) {
        VCFRotatedWindow rotation;
        vcf.getRotatedWindow({0, 3},      // Individual range
                             genomeRanges.at(i),
                             rotation,
                             lastPosition);
        Matrix partialMatrix = getMatrix(rotation);
        windowedRows = combineColumns(windowedRows, partialMatrix);
        lastPosition = rotation.posAfterLastVariant;
    }
    ASSERT_EQ(windowedRows, expectedRotation);
}