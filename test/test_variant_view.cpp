#include <gtest/gtest.h>

#include "picovcf.hpp"

#include <limits>
#include <string>

using namespace picovcf;

TEST(VariantView, Split) {
    std::string currentLine =
        "22\t1001\tv123\tA\tG,T\t.\t.\t.\tGT\t"
        "0|1\t1|1\t0|0\t1|0\t0|1";
    VCFVariantView view(currentLine);
    view.initialize(5);
    view.reset();
    ASSERT_EQ(view.getChrom(), "22");
    ASSERT_EQ(view.getPosition(), 1001);
    ASSERT_EQ(view.getID(), "v123");
    ASSERT_EQ(view.getRefAllele(), "A");
    std::vector<std::string> alts = {"G", "T"};
    ASSERT_EQ(view.getAltAlleles(), alts);
    ASSERT_TRUE(view.hasGenotypeData());
    std::vector<AlleleT> gt = {0, 1, 1, 1, 0, 0, 1, 0, 0, 1};
    ASSERT_EQ(view.getGenotypeArray(), gt);
    ASSERT_EQ(view.getMaxPloidy(), 2);
    ASSERT_EQ(view.getPhasedness(), PVCFP_PHASED);
    for (size_t i = 0; i < view.getIsPhased().size(); i++) {
        ASSERT_TRUE(view.getIsPhased()[i]);
    }
}

TEST(VariantView, Mixtures) {
    std::string currentLine =
        "arbitrary\t29854309\tv123\tA\tG,TTTTTTTTTTT\t.\t.\t.\tGT\t"
        "0/1\t1\t0|0|0\t1|0\t1|1";
    VCFVariantView view(currentLine);
    view.initialize(5);
    view.reset();
    ASSERT_EQ(view.getChrom(), "arbitrary");
    ASSERT_EQ(view.getPosition(), 29854309);
    ASSERT_EQ(view.getID(), "v123");
    ASSERT_EQ(view.getRefAllele(), "A");
    std::vector<std::string> alts = {"G", "TTTTTTTTTTT"};
    ASSERT_EQ(view.getAltAlleles(), alts);
    ASSERT_TRUE(view.hasGenotypeData());
    std::vector<AlleleT> gt = {
        0, 1, MIXED_PLOIDY,
        1, MIXED_PLOIDY, MIXED_PLOIDY,
        0, 0, 0,
        1, 0, MIXED_PLOIDY,
        1, 1, MIXED_PLOIDY,
    };
    ASSERT_EQ(view.getGenotypeArray(), gt);
    ASSERT_EQ(view.getMaxPloidy(), 3);
    ASSERT_EQ(view.getPhasedness(), PVCFP_MIXED);
    ASSERT_FALSE(view.getIsPhased().at(0));
    ASSERT_TRUE(view.getIsPhased().at(1));
    ASSERT_TRUE(view.getIsPhased().at(2));
    ASSERT_TRUE(view.getIsPhased().at(3));
    ASSERT_TRUE(view.getIsPhased().at(4));
}

TEST(VariantView, LottaAlleles) {
    std::string currentLine =
        "22\t1001\tv123\tA\tG,T,C,GG,GT,GC,TT,TC,CC,GGG,GGT,GGC,GTT,GTC,GCC,TTT,TTG,TGG,TGC\t.\t.\t.\tGT\t"
        "19|1\t10|1\t11|0\t13|.\t.|.";
    VCFVariantView view(currentLine);
    view.initialize(5);
    view.reset();
    ASSERT_EQ(view.getChrom(), "22");
    ASSERT_EQ(view.getPosition(), 1001);
    ASSERT_EQ(view.getID(), "v123");
    ASSERT_EQ(view.getRefAllele(), "A");
    std::vector<std::string> alts = {
        "G", "T", "C", "GG", "GT", "GC", "TT", "TC", "CC", "GGG", "GGT", "GGC", "GTT",
        "GTC", "GCC", "TTT", "TTG", "TGG", "TGC"};
    ASSERT_EQ(view.getAltAlleles(), alts);
    ASSERT_TRUE(view.hasGenotypeData());
    std::vector<AlleleT> gt = {
        19, 1, 10, 1, 11, 0, 13, MISSING_VALUE, MISSING_VALUE, MISSING_VALUE};
    ASSERT_EQ(view.getGenotypeArray(), gt);
    ASSERT_EQ(view.getMaxPloidy(), 2);
    ASSERT_EQ(view.getPhasedness(), PVCFP_PHASED);
    for (size_t i = 0; i < view.getIsPhased().size(); i++) {
        ASSERT_TRUE(view.getIsPhased()[i]);
    }
}