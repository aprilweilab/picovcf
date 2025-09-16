Overview of IGD
===============

An IGD file contains variant positions, identifiers (e.g., "ID" field from VCF), genotype data for a list of samples,
and identifiers for the individuals corresponding to those samples. The genotype data is referenced by *sample index*.
The sample index ``i`` ranges from ``0`` to ``N-1``, where ``N`` is the number of haplotypes. Sample indices are grouped
by individual, so given a ploidy ``P``, every ``P`` consecutive sample indices will be for the same individual.

Variants
--------

Variants in IGD are stored bi-allelic, with each IGD variant storing the sample list for a single alternate allele. So if a
site is multi-allelic with ``K`` alternate alleles, there will be ``K`` variants in the IGD file, each with the same position
and reference allele, but different alternate alleles and different sample lists.

The IGD file contains a position table, which can be scanned much faster than scanning all of the genotype data. This table
contains the position and flags for each variant, and points to the location of the genotype data. The index can be searched
linearly (:cpp:func:`picovcf::IGDData::getPosition`) or via binary search (:cpp:func:`picovcf::IGDData::lowerBoundPosition`).

Each variant (optionally) has an identifier, which can be looked up via the index of the variant in the array returned by
:cpp:func:`picovcf::IGDData::getVariantIds`.

Variant Indices
~~~~~~~~~~~~~~~

Everything variant related (genotype data, variant identifiers, reference alleles, alternate alleles) in an IGD file is
looked up by the variant index. Given ``V`` variants in the file, each variant is indexed by a number between ``0`` and ``V-1``.
The variants are ordered according to ascending base-pair position on the chromosome (IGD does not sort them - the file must
be constructed with this order). So the first polymorphic position on the chromosome (for the given dataset) is given
by ``IGDData::getPosition(0)`` and the last position is ``IGDData::getPosition(IGDData::numVariants() - 1)``.

Individuals
-----------

Given ``N`` haplotype samples and a ploidy of ``P``, there will always be ``N/P`` individuals (use :cpp:func:`picovcf::IGDData::numIndividuals`).
The (optional) identifiers for these individuals can be retrieved via :cpp:func:`picovcf::IGDData::getIndividualIds`.

Genotype Data
-------------

The genotype data is retrieved as a list of sample indexes (the ones which contain the alternate allele).

Phased
~~~~~~

For phased data, each sample index corresponds to a haplotype sample, not an individual. The :cpp:func:`picovcf::IGDData::getSamplesWithAlt`
method called on a variant index `i` will return the list of samples that have the alternate allele, which can be retrieved
via :cpp:func:`picovcf::IGDData::getAltAllele`.

Unphased
~~~~~~~~

For *unphased* data, each sample index corresponds to an individual, not a haplotype. Each variant has an additional piece
of information associated with it: the number of copies of the alternate allele that the individuals have. The number of
copies is between ``1`` and ``P`` (the ploidy), and is obtained by :cpp:func:`picovcf::IGDData::getPosition` (the three-argument
version of that function returns ``numCopies`` in the third argument). For example, with diploid individuals the homozygous
individuals are returned when ``numCopies=2`` and the heterozygous individuals are returns when ``numCopies=1``. When an
individual is homozygous in the reference allele, they will not be in any sample list (homozygous for reference is the implicit
case).  The :cpp:func:`picovcf::IGDData::getSamplesWithAlt` function is still used to retrieve the corresponding sample lists.

Contigs/Sequences
-----------------

IGD does not support contigs or sequences. The general usage is to create a single IGD file per sequence
(e.g., chromosome) rather than having many sequences in a single file. For ``igdtools``, the default is
to convert all VCF contigs into the IGD file (thus losing the "contig" information, unless you export
it as metadata - see below). A more robust approach is to either only convert VCF files with a single
contig (``igdtools --contig-require-one``) or to pick the specific sequence name you want to convert
(``igdtools --contig <contig name>``).

Metadata
--------

The IGD file itself does not contain metadata (beyond variant and individual identifiers). However, ``igdtools`` supports
exporting variant-based metadata to files that can be loaded with `numpy.loadtxt <https://numpy.org/doc/2.2/reference/generated/numpy.loadtxt.html>`_.
Matrix-based metadata (i.e., for VCF this means FORMAT fields other than GT) is not supported: if you need per-variant-per-sample metadata, then there
is probably no reason to use IGD (you need a non-sparse representation like VCF/BCF, since your metadata is non-sparse).

There are two ways to export this metadata:

1. During VCF->IGD conversion: ``igdtools in.vcf.gz -o out.igd -e all``

2. Only export metadata from a VCF: ``igdtools in.vcf.gz -e all``

The metadata is stored as a file per metadata item type. The supported fields are CHROM, QUAL, FILTER, and INFO. For INFO, each
key gets its own file.  All metadata files are a single entry (line) per variant in the resulting IGD file (i.e., "expanded" variants).

The first line of a metadata file is a comment that has information about the metadata. When loaded with ``numpy.loadtxt()``, the size of
the array is exactly :cpp:func:`picovcf::IGDData::numVariants` in length, and if you index variant ``i`` in the IGD file you can get its metadata by
looking at element ``i`` of the corresponding metadata array.

When a metadata value is not provided for a particular variant, a default value is used based on the Type field in the VCF metadata:

* Integer: ``0``
* Float: ``NaN``
* String: ``.``

Below is some example C++ code for loading metadata files. See ``examples/igd_with_meta.cpp`` for the full runnable example.
See the `pyigd documentation <https://pyigd.readthedocs.io/en/latest>`_ for a more succinct example using Python (``pyigd`` and ``numpy``).

::

  template <typename T>
  T lineToValue(const std::string& line);

  template <>
  std::string lineToValue(const std::string& line) {
      return line;
  }

  template <>
  uint64_t lineToValue(const std::string& line) {
      char* endPtr = nullptr;
      auto result = static_cast<uint64_t>(std::strtoull(line.c_str(), &endPtr, 10));
      if (endPtr == line.c_str()) {
          throw std::runtime_error("Could not parse integer");
      }
      return result;

  }

  // Function that turns a metadata file into a vector of length numVariants()
  template <typename T>
  std::vector<T> readMeta(const std::string& filename) {
      std::vector<T> result;
      std::ifstream metaTextFile(filename);
      std::string line;
      while (std::getline(metaTextFile, line)) {
          // Skip comments
          if (!line.empty() && line[0] == '#') {
              continue;
          }
          result.emplace_back(lineToValue<T>(line));
      }
      return std::move(result);
  }
