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
