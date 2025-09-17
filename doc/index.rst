picovcf Documentation
=====================

picovcf is a single-header C++ library for fast/low-memory VCF (Variant Call Format) parsing. Gzipped VCF (.vcf.gz) is also supported (optionally) if you link against zlib.

There are a lot of great tools for processing VCF files out there, but not many C++ libraries that are small (only parsing, no extra functionality) and easy to use. `picovcf` attempts to fill this niche by providing a header-only library using modern C++ (C++11) that allows clients to be selective about which parts of the VCF file get parsed.

There are two parts to the API:

1. VCF parsing.

2. Indexable Genotype Data (IGD) file format support, both reading and writing.

Example usage for both is given below.

Additionally, there is a commandline tool for converting VCF to IGD, and processing IGD files, called ``igdtools``. These tools are built when you build
picovcf, or you can also install them via ``pip install igdtools``.

.. contents::
   :depth: 2

.. toctree::
  :maxdepth: 2

  igd_overview
  igdtools
  vcf_api
  igd_api

VCF Features
------------

Not all VCF parsing features are supported. The primary focus is supporting reference-aligned
VCF data. VCF is extremely flexible, and that flexibility can make it more complex to use for
simple use cases, so ``picovcf`` tries to focus on the simpler use cases.

Features:

* Mixed ploidy and mixed phasedness are supported. See :cpp:func:`picovcf::VCFVariantView::getPhasedness` for a helper that simplifies detecting mixed phasedness.
* Variant-based metadata is supported. The main ``VCFVariantView`` accessor provides strings for the metadata (via :cpp:func:`picovcf::VCFVariantView::getMetaInfo`), and then :cpp:func:`picovcf_parse_structured_meta` is used to convert those strings into dictionaries of keys and values. See ``igdtools.cpp`` for example code.
* Non-GT per-variant, per-individual values are not supported.
* Contigs are supported. When you open a :cpp:class:`VCFFile` you choose whether you want to traverse all contigs in the file (:cpp:define:`PVCF_VCFFILE_CONTIG_ALL`), the first occurring contig (:cpp:define:`PVCF_VCFFILE_CONTIG_FIRST`), you can require the VCF to only contain a single contig (:cpp:define:`PVCF_VCFFILE_CONTIG_REQUIRE_ONE`), or you can just specify the contig name as a string and all other contigs will be ignored during traversal.

  * Note: IGD files do not support contigs. When converting from VCF to IGD you will have to decide on the same contig options above, where "traverse" becomes "convert to a single contig".

* Tabix indexes are supported on BGZF-compressed files. See :cpp:func:`picovcf::VCFFile::lowerBoundPosition` for a method that is substantially faster with indexes, and :cpp:class:`picovcf::TabixIndex` for the tabix index handling (if you want to load the index yourself).

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
      vcf.seekBeforeVariants();
      while (vcf.nextVariant()) {
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
      vcf.seekBeforeVariants();
      while (vcf.nextVariant()) {
        picovcf::VCFVariantView variant = vcf.currentVariant();
        assert(vcf.currentVariant().hasGenotypeData());
        const std::vector<AlleleT> gt = variant.getGenotypeArray();
        const size_t ploidy = variant.getMaxPloidy();
        for (size_t indiv = 0; indiv < vcf.numIndividuals(); indiv++) {
          const bool isPhased = variant.getIsPhased()[indiv];
          std::cout << (isPhased ? "phased" : "unphased") << " alleles: ";
          const size_t indivStart = indiv * ploidy;
          for (size_t j = 0; j < ploidy; j++) {
            const auto allele = gt.at(indivStart + j);
            // MIXED_PLOIDY means the current individual has a ploidy less than j
            if (allele == picovcf::MIXED_PLOIDY) {
              std::cout << "N/A, ";
            // Missing data - i.e., "." in the VCF file
            } else if (allele == picovcf::MISSING_VALUE) {
              std::cout << "?, ";
            // A regular REF/ALT value 0,1,...
            } else {
              std::cout << allele << ", ";
            }
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
 
