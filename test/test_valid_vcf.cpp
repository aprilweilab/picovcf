#include <gtest/gtest.h>

#include "picovcf.hpp"

#include <string>

using namespace picovcf;

extern const std::string getVCF42_EXAMPLE_FILE();
extern const std::string getMSPRIME_EXAMPLE_FILE();
extern const std::string getHAPLOID_DATA_EXAMPLE_FILE();
extern const std::string getMIXED_DATA_EXAMPLE_FILE();
extern const std::string getCONTIG_EXAMPLE_FILE();

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
        ASSERT_TRUE(vcf.nextVariant());
    }
    ASSERT_FALSE(vcf.nextVariant());

    vcf.seekBeforeVariants();
    ASSERT_TRUE(vcf.nextVariant());
    VCFVariantView& variant1 = vcf.currentVariant();
    ASSERT_EQ(variant1.getChrom(), "20");
    ASSERT_EQ(variant1.getPosition(), 14370);
    ASSERT_EQ(variant1.getID(), "rs6054257");
    ASSERT_EQ(variant1.getRefAllele(), "G");
    ASSERT_EQ(variant1.getAltAlleles(), std::vector<std::string>{"A"});
    ASSERT_TRUE(variant1.hasGenotypeData());
    ASSERT_EQ(vcf.numIndividuals(), 3);

    // Check all individuals of the first variant
    VariantT allele1 = MISSING_VALUE;
    VariantT allele2 = MISSING_VALUE;
    bool isPhased;
    const auto& gtArray = variant1.getGenotypeArray();
    const size_t ploidy = variant1.getMaxPloidy();
    ASSERT_EQ(variant1.getPhasedness(), PVCFP_MIXED);
    ASSERT_EQ(ploidy, 2);
    // The first individual.
    ASSERT_TRUE(variant1.getIsPhased().at(0));
    ASSERT_EQ(gtArray.at(0), 0);
    ASSERT_EQ(gtArray.at(1), 0);
    // Second individual.
    ASSERT_TRUE(variant1.getIsPhased().at(1));
    ASSERT_EQ(gtArray.at(2), 1);
    ASSERT_EQ(gtArray.at(3), 0);
    // Third individual.
    ASSERT_FALSE(variant1.getIsPhased().at(2));
    ASSERT_EQ(gtArray.at(4), 1);
    ASSERT_EQ(gtArray.at(5), 1);

    // Check a few things on the 4th variant
    ASSERT_TRUE(vcf.nextVariant());
    ASSERT_TRUE(vcf.nextVariant());
    ASSERT_TRUE(vcf.nextVariant());
    VCFVariantView& variant4 = vcf.currentVariant();
    ASSERT_EQ(variant4.getChrom(), "20");
    ASSERT_EQ(variant4.getPosition(), 1230237);
    ASSERT_EQ(variant4.getID(), ".");
    ASSERT_EQ(variant4.getRefAllele(), "T");
    ASSERT_EQ(variant4.getAltAlleles(), std::vector<std::string>());
    ASSERT_TRUE(variant4.hasGenotypeData());
}

static bool hasSample(const IGDSampleList& samples, SampleT id) {
    for (auto sampleId : samples) {
        if (id == sampleId) {
            return true;
        }
    }
    return false;
}

// Equality between VCF and IGD on the same data.
TEST(ValidVCF, Indexable) {
    constexpr size_t EXPECT_VARIANTS = 4;
    constexpr size_t EXPECT_MISS = 3;
    constexpr size_t EXPECT_INDIVIDUALS = 10000;

    // Verify the VCF file first.
    VCFFile vcf(getMSPRIME_EXAMPLE_FILE());
    ASSERT_EQ(vcf.numVariants(), EXPECT_VARIANTS);
    ASSERT_EQ(vcf.numIndividuals(), EXPECT_INDIVIDUALS);

    // Now convert to Indexable Genotype Data
    std::string testFilename = "tmp.indexable.igd";
    vcfToIGD(getMSPRIME_EXAMPLE_FILE(), testFilename);
    IGDData indexableData(testFilename);
    ASSERT_EQ(indexableData.numVariants(), EXPECT_VARIANTS+EXPECT_MISS);
    ASSERT_EQ(indexableData.numIndividuals(), EXPECT_INDIVIDUALS);
    ASSERT_EQ(indexableData.numSamples(), 2*EXPECT_INDIVIDUALS);
    bool _ignore = false;
    ASSERT_EQ(indexableData.getPosition(0, _ignore), 55829);

    size_t index = 0;
    vcf.seekBeforeVariants();
    while (vcf.nextVariant()) {
        VCFVariantView& variant = vcf.currentVariant();
        ASSERT_FALSE(variant.getAltAlleles().size() > 1);
        bool isMissing = false;
        const auto position = indexableData.getPosition(index, isMissing);
        if (isMissing) {
            continue;
        }
        ASSERT_EQ(variant.getPosition(), position);
        auto sampleSet = indexableData.getSamplesWithAlt(index);
        const auto& gtArray = variant.getGenotypeArray();
        const size_t ploidy = variant.getMaxPloidy();
        ASSERT_LE(ploidy, 2);
        for (size_t indiv = 0; indiv < vcf.numIndividuals(); indiv++) {
            ASSERT_TRUE(variant.getIsPhased().at(indiv));
            const size_t baseIndex = indiv*ploidy;
            for (size_t j = 0; j < ploidy; j++) {
                const SampleT sampleIndex = baseIndex + j;
                if (gtArray.at(sampleIndex) > 0 && gtArray.at(sampleIndex) < MISSING_VALUE) {
                    ASSERT_TRUE(hasSample(sampleSet, sampleIndex));
                }
            }
        }
        index++;
    }
}

TEST(ValidVCF, Haploid) {
    const size_t EXPECT_VARIANTS = 4;
    VCFFile vcf(getHAPLOID_DATA_EXAMPLE_FILE());

    ASSERT_EQ(vcf.numVariants(), EXPECT_VARIANTS);
    ASSERT_EQ(vcf.numIndividuals(), 10000);

    vcf.seekBeforeVariants();
    for (size_t i = 0; i < EXPECT_VARIANTS; i++) {
        ASSERT_TRUE(vcf.nextVariant());
        VCFVariantView variant = vcf.currentVariant();
        const auto& gtArray = variant.getGenotypeArray();
        const size_t ploidy = variant.getMaxPloidy();
        ASSERT_EQ(ploidy, 1);
        ASSERT_EQ(variant.getPhasedness(), PVCFP_PHASED);
        for (size_t indiv = 0; indiv < vcf.numIndividuals(); indiv++) {
            ASSERT_TRUE(variant.getIsPhased().at(indiv));
        }
    }
    ASSERT_FALSE(vcf.nextVariant());
}

TEST(ValidVCF, MixedPloidy) {
    const size_t EXPECT_VARIANTS = 4;
    VCFFile vcf(getMIXED_DATA_EXAMPLE_FILE());

    ASSERT_EQ(vcf.numVariants(), EXPECT_VARIANTS);
    ASSERT_EQ(vcf.numIndividuals(), 10000);

    vcf.seekBeforeVariants();
    for (size_t i = 0; i < EXPECT_VARIANTS; i++) {
        ASSERT_TRUE(vcf.nextVariant());
        VCFVariantView variant = vcf.currentVariant();
        const auto& gtArray = variant.getGenotypeArray();
        const size_t ploidy = variant.getMaxPloidy();
        ASSERT_EQ((double)gtArray.size() / (double)ploidy, (double)vcf.numIndividuals());
        if (i < 2) {
            ASSERT_EQ(ploidy, 1);
        } else {
            ASSERT_EQ(ploidy, 2);
        }
        for (size_t indiv = 0; indiv < vcf.numIndividuals(); indiv++) {
            ASSERT_TRUE(variant.getIsPhased().at(indiv));
            const size_t baseIndex = indiv*ploidy;
            ASSERT_NE(gtArray.at(baseIndex), NOT_DIPLOID);
            if (i >= 2) {
                ASSERT_NE(gtArray.at(baseIndex+1), NOT_DIPLOID);
            }
        }
    }
    ASSERT_FALSE(vcf.nextVariant());
}

std::vector<size_t> getAllPositions(VCFFile& vcf) {
    std::vector<size_t> result;
    vcf.seekBeforeVariants();
    while (vcf.nextVariant()) {
        VCFVariantView& variant = vcf.currentVariant();
        result.push_back(variant.getPosition());
    }
    return std::move(result);
}

TEST(ValidVCF, Contigs) {
    const std::string inputFile = getCONTIG_EXAMPLE_FILE();
    // All contigs (default)
    {
        const size_t EXPECT_VARIANTS = 6;
        VCFFile vcf(inputFile);

        ASSERT_EQ(vcf.numVariants(), EXPECT_VARIANTS);
        ASSERT_EQ(vcf.numIndividuals(), 3);
        auto positions = getAllPositions(vcf);
        const std::vector<size_t> expected = {4000, 4500, 7777, 14000, 14500, 17777};
        ASSERT_EQ(positions, expected);
    }

    // Only one contig -- should fail
    {
        VCFFile vcf(inputFile, PVCF_VCFFILE_CONTIG_REQUIRE_ONE);
        ASSERT_THROW(vcf.numVariants(), ApiMisuse);
    }

    // Non-existent contig -- should fail
    {
        ASSERT_THROW(
            VCFFile vcf(inputFile, "random_string"),
            ApiMisuse);
    }

    // contig "20" by name and by being the first
    {
        const size_t EXPECT_VARIANTS = 3;
        VCFFile vcf(inputFile, "20");
        VCFFile vcfFirst(inputFile, PVCF_VCFFILE_CONTIG_FIRST);

        ASSERT_EQ(vcf.numVariants(), EXPECT_VARIANTS);
        ASSERT_EQ(vcfFirst.numVariants(), EXPECT_VARIANTS);
        ASSERT_EQ(vcf.numIndividuals(), 3);
        ASSERT_EQ(vcfFirst.numIndividuals(), 3);
        const auto positions = getAllPositions(vcf);
        const auto firstPos = getAllPositions(vcfFirst);
        const std::vector<size_t> expected = {4000, 4500, 7777};
        ASSERT_EQ(positions, expected);
        ASSERT_EQ(firstPos, expected);
    }

    // contig "blargh" by name
    {
        const size_t EXPECT_VARIANTS = 3;
        VCFFile vcf(inputFile, "blargh");

        ASSERT_EQ(vcf.numVariants(), EXPECT_VARIANTS);
        ASSERT_EQ(vcf.numIndividuals(), 3);
        const auto positions = getAllPositions(vcf);
        const std::vector<size_t> expected = {14000, 14500, 17777};
        ASSERT_EQ(positions, expected);
    }
}
