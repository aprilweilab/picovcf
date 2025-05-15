# igdtools

`igdtools` can convert from .vcf(.gz) to IGD, and once you have an IGD file it can perform
various operations such as filtering, computing basic statistics, and generally transforming
IGD files.

Run `igdtools --help` for more information on commands.

For more general reading and modification of IGD files, see [pyigd](https://github.com/aprilweilab/pyigd).

## Installation

`igdtools` is a C++ binary, with a small Python wrapper to make installation easier. You can install
via:

```
pip install igdtools
```

which will install prebuilt binaries for most Linux systems, and install via a source distribution
for other systems (such as MacOS). The source distribution requires CMake 3.10 or newer, the
zlib development headers, and a version of clang or GCC that supports C++11.

## Usage Examples

### Convert .vcf(.gz) to IGD

Conversion will copy the variant identifiers ("ID" column in VCF) and individual identifiers (the
sample column names in VCF) to the IGD file, unless `--no-var-ids` and `--no-indiv-ids` flags
are specified (respectively).

```
igdtools input.vcf.gz -o output.igd
```

### Convert .vcf(.gz) to IGD and export metadata

`igdtools` can export metadata fields as simple text files, each of which can be
loaded by [numpy.loadtxt()](https://numpy.org/doc/2.2/reference/generated/numpy.loadtxt.html).
You can use ``--export-metadata`` to export this metadata during conversion:

```
igdtools input.vcf.gz -o output.igd --export-metadata qual,filter,info
```

The list of metadata types you can export are `qual,filter,chrom,info`. You can also
specify `all` to just export all of them without listing them out. By default,
no metadata is exported during VCF to IGD conversion.

### Just export .vcf(.gz) metadata

If you already have an IGD file, and want to go back and export the metadata
from the corresponding VCF, you can just do the export:

```
igdtools input.vcf.gz --export-metadata qual,filter,info
```

Note that the naming will differ. When just exporting the metadata, the naming
follows `input.meta.*.txt`, whereas the previous example with convering to IGD
_and_ exporting metadata would have named based on `output.meta.*.txt`.

### IGD file header info

To examine the header information in the IGD file:

```
igdtools -i test.igd
```

which will print out something like
```
  Variants: 329556
  Individuals: 1000
  Ploidy: 2
  Phased?: true
  Source: true_data/simulation-source-1000-100mb.vcf
  Genome range: 115-99999629
  Has individual IDs? Yes
  Has variant IDs? No
```

Similarly, some simple statistics can be emitted by specifying `-s`, which causes
the entire IGD file to be scanned (so will be slower than `-i`).

### Copy variants within range

Create a copy of an IGD file, but only keep variants in a particular base-pair
range. Here we show range `[50000, 150000]` (50KB to 150KB, inclusive):

```
igdtools input.igd -o output.igd -r 50000-150000
```

### Copy variants with frequency

Create a copy of an IGD file, but only keep variants that have a particular
allele frequency range. Here we show range `[0.1, 0.4)` (inclusive, exclusive).

```
igdtools input.igd -o output.igd -f 0.1-0.4
```

### Copy to unphased data

Sometimes it is useful to perform "unphased" calculations. For example, when
computing runs of homozygosity (ROH) it is easier to work with unphased diploid
data that tracks the number of copies each individual has of an allele (0, 1, or 2).
Create a copy of an IGD file, but store it unphased:

```
igdtools test.igd -o test.unphased.igd --force-unphased
```