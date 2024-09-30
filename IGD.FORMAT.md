# IGD Binary Format

IGD has a fixed-size header that contains offsets to other parts of the file. The header
starts at the 0th byte of the file, and is 128 bytes in size:

| Byte Offset | Byte Size | Name | Description |
| ---- | ---- | ---- | ---- |
| 0 | 8 | magic | Magic number that identifies the file as an IGD: 0x3a0c6fd7945a3481 |
| 8 | 8 | version | The version of the file format. Only incremented on breaking changes |
| 16 | 4 | ploidy | Ploidy of each individual, >= 1 |
| 20 | 4 | sparseThreshold | Threshold that determines whether a sample list is stored as a bitvector or sparse list. |
| 24 | 8 | numVariants | The number of total variants in the file |
| 32 | 4 | numIndividuals | The number of total individuals in the file |
| 36 | 4 | reserved | Reserved for future field. MUST be set to 0 when an IGD file is created |
| 40 | 8 | flags | Bitwise flags; the only one right now is 0x1, which if set means the data is phased |
| 48 | 8 | filePosIndex | The file position of the index describing the genomic and file positions of all variants. |
| 56 | 8 | filePosVariants | The file position of the first variant data entry - this data does not contain genotypes, just information about variants |
| 64 | 8 | filePosIndividualIds | The file position of the identifiers (labels) for individuals, set to 0 if there are no identifiers |
| 72 | 8 | filePosVariantIds | The file position of the identifiers (labels) for variants, set to 0 if there are no identifiers |
| 80 | 48 | reserved | Reserved data for future fields. MUST be set to 0 when an IGD file is created |

All variable-sized strings in the file format are encoded as an 4-byte length `L`, followed by `L` bytes containing the string contents.
All sample-related quantities are 32-bit (allowing a maximum of ~4 billion samples, ~2 billion individuals). All variant-related quantities are 64-bit.
All indexes are 0-based, so `0` represents the first sample, first variant, etc.

Immediately after the header are two variable-sized string fields:
1. "Source": where the IGD file data came from (often information about data conversion).
2. "Description": just a general field to hold information about the dataset.

## Sample Rows (REQUIRED)

Immediately following the header and string fields is the sample row data. This data should be accessed only via the corresponding `IndexEntry`, which describes the layout.

Note that we only support one alternate per variant, which means that converting from a file format like VCF requires "expanding" the variants. A single VCF variant with `N` alternate alleles will become `N` IGD variants.

The header, string firleds, and sample rows are the only required-to-be-contiguous parts of the IGD file. The remaining sections can be anywhere in the file (_after_ the above), and use the header fields to indicate where they are.
Each row of sample data corresonds to the list of samples that have the alternate allele of the given (bi-allelic) variant. When the row is flagged as a "missing data row", the sample list corresponds to samples that do not have an allele for the given site, in the dataset. Each row of sample data can be encoded as either a "SizedList" or "BitVector" -- see below.

### SizedList Samples

The "SizedList" is a sparse representation of the samples as a list of indexes. If the list contains value `i`, then the `ith` sample has the alternate allele represented by the particular row.

The layout is:
  * A 4-byte unsigned integer: the number of samples `k` (NOT individuals) that are in the list.
  * `k` consecutive 4-byte unsigned integers: the index that corresponds to the sample in the list: this value is between `0...((ploidy*numIndividuals)-1)`.

### BitVector Samples

The non-sparse representation of the samples as a bitvector. A `1` indicates the alternate allele, a `0` represents the reference allele. There are `ploidy*numIndividuals` bits, rounded up to the next byte (i.e., `ceil( (ploidy*numIndividuals) / 8)`). The end of the row can by found by `startingPosition + (((numIndividuals * ploidy) + 7) / 8)`. Bit `startingPosition + i` refers to the `ith` sample.

## Variant Information (REQUIRED)

The string values for alleles are located in the file at `filePosVariants` in the header. Has `numVariants` entries, each of which is a pair of string:
* Variable-sized string: Reference allele value.
* Variable-sized string: Alternate allele value.

Note that we only support one alternate per variant, which means that converting from a file format like VCF requires "expanding" the variants. A single VCF variant with `N` alternate alleles will become `N` IGD variants.

## Individual Identifiers (OPTIONAL)

Identifiers for each individual are located in the file at `filePosIndividualIds` if it is not `0`. This data has:
* An 8-byte value indicating how many labels there are.
* A variable-sized string for each label.

## Variant Identifiers (OPTIONAL)

Identifiers for each variant are located in the file at `filePosVariantIds` if it is not `0`. This data has:
* An 8-byte value indicating how many labels there are.
* A variable-sized string for each label.

## The Index (REQUIRED)

At position identified by `filePosIndex` is the file index. There are `numVariants` consecutive `IndexEntry`s following this layout:
| Byte Offset | Byte Size | Name | Description |
| ---- | ---- | ---- | ---- |
| 0 | 1 | flags | Flags describing properties of this variant. Can be 0 (nothing) or a bitwise combination of: `SPARSE=0x1`, `IS_MISSING=0x2` |
| 1 | 7 | bpPosition | The position on the genome in base pairs for this variant. |
| 8 | 8 | filePosDataRow | The file position for the sample row. This can be sparse (see `SizedList`) or a bit vector (see `BitVector`), depending on if flag `SPARSE` is set. If `IS_MISSING` flag is set, this represents missing data from this variant. |
