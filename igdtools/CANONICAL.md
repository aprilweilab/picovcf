# Canonical IGD

The IGD format is a fairly flexible representation of sparse matrices. There are very few restrictions
on the relationships between rows (variants) in the dataset. This flexibility is useful for a file format,
but can introduce complexities for downstream processing because such tools/methods need to consider many
different scenarios.

This document describes the "canonical" IGD format, which imposes some further restrictions on IGD files
that simplify downstream processing. Canonical IGD is not a change to the format, but can be thought of
as a "check" on an IGD file that either passes (it is canonical) or fails (it is not canonical).

## Standard IGD Constraints

1. Variants (rows) are in position-ascending sorted order
2. No variant (row) contains more than `N` samples, where `N` is the number of individuals (unphased data) or haplotype samples (phased data).
3. The minimum position value is `0` and the maximum is `(2^64) - 3`.

## Canonical Constraints

### Non-overlapping Site Sample Sets

Consider a site at position `P`. There can be an arbitrary set of variants at that site, `V`. The following
are the constraints on this set of variants `V`:
1. There is a single REF allele at the site (for all SNPs _and_ non-SNPs like indels, etc.)
2. There is a single variant (row) for each alternate allele
3. Let `S_i` be the list of sample identifiers for the `i`th variant in `V`. Then the intersection of `S_i` and `S_j` must be empty for all `i` not equal `j` (and `1 <= j <= |V|`).
4. There is a single variant (row) for the missing data at the site. The missing data sample list also follows the empty intersection rule above.

### Non-SNPs (indels, etc)

There are no additional constraints on non-SNPs.

### Unphased Data

There are two main differences of unphased data from phased data, in IGD:
1. The sample lists are made up of _individual_ indexes, instead of _haplotype_ indexes.
2. Each variant is associated with a value `NumCopies`. For a ploidy of `P >= 1`, `0 <= NumCopies <= P`.
 
The constraints on unphased data are the same as for phased data, except:
1. `NumCopies=0` is reserved for use with missing data only.
2. There may be multiple variants (rows) for each alterante allele, but _only if they have different values for `NumCopies`_.

## Creating a canonical IGD

Given a non-canonical IGD, how do we convert it to a canonical IGD? In some cases, the user should apply their own
QA-related rules and investigate any remaining violations of canonical rules, and resolve them manually.

Another option is to just drop all variants at any site that violates canonical IGD rules. As long as the number of violations is small (proportionally), this may be the best approach. However,
filtering by removing sites can affect analyses that rely on segregating sites between samples (which includes many
population genetic analyses).

If we want an automatic way to create a canonical IGD without removing all violating sites, we need a deterministic
set of rules to convert violating sites into canonical IGD compliant sites. We also want an audit trail of all
changes made.

Below are some rules for conversion if the user does not want to drop the entire site(s).

### Rule: single REF allele at the site

When there are multiple REF alleles, the alleles are sorted according to [shortlex order](https://en.wikipedia.org/wiki/Shortlex_order), so SNPs are ordered before non-SNPs, and then the first allele is chosen for a reference. The variants associated with non-chosen reference alleles are dropped.

### Rule: single variant (row) for each alternate allele

The sample sets of the two matching alternate alleles are merged.

### Rule: empty intersection of sample sets

Consider every pair of sample sets `S_i`, `S_j` that overlap.
* If one variant is missing data, let it be `i`, and then set `S_i = S_i \ S_j` (i.e., remove the samples from the missing data)
* Otherwise, create a new variant with the union of the sample sets, with alternate allele `<A1>_OR_<A2>` where `<A1>`, `<A2>` are the original alleles.