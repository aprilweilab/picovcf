#include <gtest/gtest.h>

#include "picovcf.hpp"

#include <string>

using namespace picovcf;

extern const std::string getVCF42_EXAMPLE_FILE();
extern const std::string getMSPRIME_EXAMPLE_FILE();

TEST(ValidVCF, SpecExample) {
    const size_t EXPECT_VARIANTS = 6;
    VCFFile vcf(getVCF42_EXAMPLE_FILE());

    std::string phasing = vcf.getMetaInfo("phasing");
    ASSERT_EQ(phasing, "partial");
    std::string source = vcf.getMetaInfo("source");
    ASSERT_EQ(source, "myImputationProgramV3.1");

    ASSERT_EQ(vcf.numVariants(), EXPECT_VARIANTS);
    ASSERT_EQ(vcf.numIndividuals(), 3);
    ASSERT_EQ(vcf.getIndividualLabels(), std::vector<std::string>({"NA00001", "NA00002", "NA00003"}));

    vcf.seekBeforeVariants();
    for (size_t i = 0; i < EXPECT_VARIANTS; i++) {
        ASSERT_TRUE(vcf.hasNextVariant());
        vcf.nextVariant();
        VCFVariantView variant = vcf.currentVariant();
        IndividualIteratorGT iterator = variant.getIndividualIterator();
        while (iterator.hasNext()) {
            IndexT allele1 = 0;
            IndexT allele2 = 0;
            iterator.getAlleles(allele1, allele2);
        }
    }
    ASSERT_FALSE(vcf.hasNextVariant());

    vcf.seekBeforeVariants();
    ASSERT_TRUE(vcf.hasNextVariant());
    vcf.nextVariant();
    VCFVariantView& variant1 = vcf.currentVariant();
    ASSERT_EQ(variant1.getChrom(), "20");
    ASSERT_EQ(variant1.getPosition(), 14370);
    ASSERT_EQ(variant1.getID(), "rs6054257");
    ASSERT_EQ(variant1.getRefAllele(), "G");
    ASSERT_EQ(variant1.getAltAlleles(), std::vector<std::string>{"A"});
    ASSERT_TRUE(variant1.hasGenotypeData());

    // Check all individuals of the first variant
    size_t allele1 = MISSING_VALUE;
    size_t allele2 = MISSING_VALUE;
    bool isPhased;
    IndividualIteratorGT iterator = variant1.getIndividualIterator();
    // Do the first individual twice in a row.
    ASSERT_TRUE(iterator.hasNext());
    isPhased = iterator.getAlleles(allele1, allele2, /*moveNext=*/false);
    ASSERT_TRUE(isPhased);
    ASSERT_EQ(allele1, 0);
    ASSERT_EQ(allele2, 0);
    isPhased = iterator.getAlleles(allele1, allele2);
    ASSERT_TRUE(isPhased);
    ASSERT_EQ(allele1, 0);
    ASSERT_EQ(allele2, 0);

    ASSERT_TRUE(iterator.hasNext());
    isPhased = iterator.getAlleles(allele1, allele2);
    ASSERT_TRUE(isPhased);
    ASSERT_EQ(allele1, 1);
    ASSERT_EQ(allele2, 0);

    ASSERT_TRUE(iterator.hasNext());
    isPhased = iterator.getAlleles(allele1, allele2);
    ASSERT_FALSE(isPhased);
    ASSERT_EQ(allele1, 1);
    ASSERT_EQ(allele2, 1);

    ASSERT_FALSE(iterator.hasNext());

    // Check a few things on the 4th variant
    ASSERT_TRUE(vcf.hasNextVariant());
    vcf.nextVariant();
    ASSERT_TRUE(vcf.hasNextVariant());
    vcf.nextVariant();
    ASSERT_TRUE(vcf.hasNextVariant());
    vcf.nextVariant();
    VCFVariantView& variant4 = vcf.currentVariant();
    ASSERT_EQ(variant4.getChrom(), "20");
    ASSERT_EQ(variant4.getPosition(), 1230237);
    ASSERT_EQ(variant4.getID(), ".");
    ASSERT_EQ(variant4.getRefAllele(), "T");
    ASSERT_EQ(variant4.getAltAlleles(), std::vector<std::string>());
    ASSERT_TRUE(variant4.hasGenotypeData());
}

using SampleList = std::vector<IndexT>;

static bool hasSample(const SampleList& samples, IndexT id) {
    for (auto sampleId : samples) {
        if (id == sampleId) {
            return true;
        }
    }
    return false;
}

TEST(ValidVCF, Indexable) {
    constexpr size_t EXPECT_VARIANTS = 4;
    constexpr size_t EXPECT_INDIVIDUALS = 10000;

    // Verify the VCF file first.
    VCFFile vcf(getMSPRIME_EXAMPLE_FILE());
    ASSERT_EQ(vcf.numVariants(), EXPECT_VARIANTS);
    ASSERT_EQ(vcf.numIndividuals(), EXPECT_INDIVIDUALS);

    // Now convert to Indexable Genotype Data
    std::string testFilename = "tmp.indexable.igd";
    vcfToIGD(getMSPRIME_EXAMPLE_FILE(), testFilename);
    IGDData indexableData(testFilename);
    ASSERT_EQ(indexableData.numVariants(), EXPECT_VARIANTS);
    ASSERT_EQ(indexableData.numIndividuals(), EXPECT_INDIVIDUALS);
    ASSERT_EQ(indexableData.numSamples(), 2*EXPECT_INDIVIDUALS);
    ASSERT_EQ(indexableData.getPosition(0), 55829);

    size_t index = 0;
    vcf.seekBeforeVariants();
    while (vcf.hasNextVariant()) {
        vcf.nextVariant();
        VCFVariantView& variant = vcf.currentVariant();
        ASSERT_FALSE(variant.getAltAlleles().size() > 1);
        ASSERT_EQ(variant.getPosition(), indexableData.getPosition(index));
        auto sampleSet = indexableData.getSamplesWithAlt(index);
        IndividualIteratorGT individualIt = variant.getIndividualIterator();
        size_t sampleIndex = 0;
        while (individualIt.hasNext()) {
            IndexT allele1 = 255;
            IndexT allele2 = 255;
            individualIt.getAlleles(allele1, allele2);
            if (allele1 == 1) {
                ASSERT_TRUE(hasSample(sampleSet, sampleIndex));
            }
            sampleIndex++;
            if (allele2 != NOT_DIPLOID) {
                if (allele2 == 1) {
                    ASSERT_TRUE(hasSample(sampleSet, sampleIndex));
                }
                sampleIndex++;
            }
        }
        index++;
    }
}
