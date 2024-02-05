# IGD Binary Format

IGD has a fixed-size header that contains offsets to other parts of the file. The header
starts at the 0th byte of the file, and is 128 bytes in size:

| Byte Offset | Byte Size | Name | Description |
| ---- | ---- | ---- | ---- |
| 0 | 8 | magic | Magic number that identifies the file as an IGD: 0x3a0c6fd7945a3481 |
| 8 | 8 | version | The version of the file format. Only incremented on breaking changes |
| 16 | 8 | ploidy | Ploidy of each individual, >= 1 |
| 24 | 8 | numVariants | The number of total variants in the file |
| 32 | 8 | numIndividuals | The number of total individuals in the file |
| 40 | 8 | flags | Bitwise flags; the only one right now is 0x1, which if set means the data is phased |
| 48 | 8 | filePosVariants | The file position of the first variant data entry - this data does not contain genotypes, just information about variants |
| 56 | 8 | filePosMissingData | The file position of the missing data map, which is an optionally-loadable part of the file |
| 64 | 8 | filePosIndividualIds | The file position of the identifiers (labels) for individuals, left at 0 if there are no identifiers |
| 72 | 56 | reserved | Reserved data for future fields. MUST be set to 0 when an IGD file is created |

All variable-sized strings in the file format are encoded as an 8-byte length `L`, followed by `L` bytes containing the string contents.

Immediately after the header are two variable-sized string fields:
1. "Source": where the IGD file data came from (often information about data conversion).
2. "Description": just a general field to hold information about the dataset.

Next is the genotype data:
* There are `numVariants` consecutive entries.
* Each entry is an 8-byte genome position, followed by `ploidy*numIndividuals` bits, rounded up to the next byte.
* Given a variant index (0-based) `i`, we can then find the genotype for that variant by starting at the 0th one and adding `i * (((numIndividuals * ploidy) + 7) / 8)` to it.

The above are the only required-to-be-consecutive parts of the IGD file. There are three remaining sections that can be anywhere in the file (_after_ the above), and use the header fields to indicate where they are.

## Variant Information

At position identified by `filePosVariants` in the header. Has `numVariants` entries, each of which is:
* Variable-sized string: Reference allele value.
* Variable-sized string: Alternate allele value.

Note that we only support one alternate per variant, which means that converting from a file format like VCF requires "expanding" the variants. A single VCF variant with `N` alternate alleles will become `N` IGD variants.

## Individual Identifiers

Identifiers for each individual are at `filePosIndividualIds` if it is not `0`. This data has:
* An 8-byte value indicating how many labels there are.
* A variable-sized string for each label.

## Missing Data

The expectation is that missing data is fairly sparse, so we store it separate from the other genotype data. The missing data is stored as:
* An 8-byte value indicating how many entries there are.
* For each entry:
  * The 8-byte genomic position that is missing the data.
  * The 8-byte number of samples (NOT individuals) that are missing the data.
  * An 8-byte index that corresponds to the sample that has missing data: this value is between `0...(ploidy*numIndividuals)`.
