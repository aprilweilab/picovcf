picovcf Documentation
=====================

picovcf is a single-header C++ library for fast/low-memory VCF (Variant Call Format) parsing. Gzipped VCF (.vcf.gz) is also supported (optionally) if you link against zlib.

There are a lot of great tools for processing VCF files out there, but not many C++ libraries that are small (only parsing, no extra functionality) and easy to use. `picovcf` attempts to fill this niche by providing a header-only library using modern C++ (C++11) that allows clients to be selective about which parts of the VCF file get parsed.


.. contents::
   :depth: 2

.. toctree::
  :maxdepth: 2

  vcf_api
  igd_api

Simple Usage Examples
---------------------

VCF Example: Iterate variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  #include <iostream>
  #include "picovcf.hpp"

  int main(int argc, char *argv[]) {
    if (argc >= 2) {
      picovcf::VCFFile vcf(argv[1]);
      vcf.firstVariant();
      while (vcf.hasNextVariant()) {
        vcf.nextVariant();
        picovcf::VCFVariantView variant = vcf.currentVariant();
        std::cout << "Variant at position: " << variant.getPosition() << std::endl;
      }
    }
    return 0;
  }


VCF Example: Iterate genotypes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  #include <iostream>
  #include <cassert>
  #include "picovcf.hpp"

  int main(int argc, char *argv[]) {
    if (argc >= 2) {
      picovcf::VCFFile vcf(argv[1]);
      vcf.firstVariant();
      while (vcf.hasNextVariant()) {
        vcf.nextVariant();
        picovcf::VCFVariantView variant = vcf.currentVariant();
        assert(vcf.currentVariant().hasGenotypeData());
        picovcf::IndividualIteratorGT iterator = variant.getIndividualIterator();
        while (iterator.hasNext()) {
          picovcf::IndexT allele1 = 0;
          picovcf::IndexT allele2 = 0;
          bool isPhased = iterator.getAlleles(allele1, allele2);
          std::cout << (isPhased ? "phased" : "unphased") << " alleles: " << allele1;
          if (allele2 != picovcf::NOT_DIPLOID) {
            std::cout << ", " << allele2;
          }
          std::cout << std::endl;
        }
      }
    }
    return 0;
  }


IGD Example: Iterate variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  #include <iostream>
  #include "picovcf.hpp"

  int main(int argc, char *argv[]) {
    if (argc >= 2) {
      picovcf::IGDData igd(argv[1]);
      for (size_t i = 0; i < igd.numVariants(); i++) {
        std::cout << "Variant at position: " << igd.getPosition(i) << std::endl;
      }
    }
    return 0;
  }
 
IGD Example: Iterate genotypes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  #include <iostream>
  #include "picovcf.hpp"

  int main(int argc, char *argv[]) {
    if (argc >= 2) {
      picovcf::IGDData igd(argv[1]);
      for (size_t i = 0; i < igd.numVariants(); i++) {
        std::cout << "Samples with alternate allele \"" << igd.getAltAllele(i) << "\": ";
        auto sampleList = igd.getSamplesWithAlt(i);
        for (picovcf::IndexT sampleIndex : sampleList) {
          std::cout << sampleIndex << ", ";
        }
        std::cout << std::endl;
      }
    }
    return 0;
  }
 
