![](https://github.com/aprilweilab/picovcf/actions/workflows/cmake-multi-platform.yml/badge.svg)
![](https://readthedocs.org/projects/picovcf/badge/?version=latest)

# picovcf

Single-header C++ library for fast/low-memory VCF (Variant Call Format) parsing. Gzipped VCF (.vcf.gz) is optionally supported.

There are a lot of great tools for processing VCF files out there, but not many C++ libraries that are small (only parsing, no extra functionality) and easy to use. `picovcf` attempts to fill this niche by providing a header-only library using modern C++ (C++11) that allows clients to be selective about which parts of the VCF file get parsed.

## Using the library

Either copy the latest header file (`picovcf.hpp`) into your project directly, or make use of something like git submodules to include https://github.com/aprilweilab/picovcf.

See the [vcfpp.cpp](https://github.com/aprilweilab/picovcf/blob/main/examples/vcfpp.cpp) for an example of how to use the APIs. Read [the docs](https://picovcf.readthedocs.io/en/latest/) for an overview of the API.

When building code that uses `picovcf.hpp`, define `VCF_GZ_SUPPORT=1` (`-DVCF_GZ_SUPPORT=1` on most compiler command lines) to enable zlib support for compressed VCF files.

## Build and run the tests

picovcf does not need to be built to be used, since it is a single header that gets built as part of your project. However, if you want to build and run the unit tests, do the following:

```
cd picovcf
mkdir build && cd build
cmake ..
make
EXAMPLE_VCFS=../test/example_vcfs/ ./picovcf_test
```

There is a Dockerfile that encodes all the build steps and dependencies, including documentation build.

## Build the documentation

Requires Python packages `sphinx`, `sphinx-rtd-theme`, `breathe`. Requires Doxygen.

From the same `build/` directory as above:
```
DOC_BUILD_DIR=$PWD sphinx-build -c ../doc/ -b html -Dbreathe_projects.picovcf=$PWD/doc/xml ../doc/ $PWD/doc/sphinx/
```

## Indexable Genotype Data (IGD)

`picovcf` also defines an extremely simple binary file format that can be used for fast access to genotype data. Most other genotype data formats are not indexable directly: that is, you cannot jump directly to the 1 millionth variant without first scanning all the previous (almost million) variants. IGD has the following properties:
* Indexable. You can use math to figure out where the `i`th variant will be in the file.
* Uncompressed. No need to link in compression libraries.
* Simple format: all variants are expanded into binary variants. So if a Variant has `N` alternate alleles, then IGD will store that as `N` rows containing `0` (reference allele) or `1` (alternate allele).
* Reasonably small. It is slightly larger than more complex file formats (like BGEN), but usually within a factor of `3`-`4` in filesize.

For example, the following are from chromosome 22 of a real dataset:
* `.vcf`: 11GB
* `.vcf.gz`: 203MB
* `.bgen`: 256MB
* `.igd`: 691MB

Converting the `.vcf.gz` to `.bgen` (via qctool) took 23 minutes, but converting to `.igd` only took 3 minutes. Furthermore, iteratively accessing all the variants (and genotype data) in the `.igd` file was ~`5x` faster than accessing the same data in the `.vcf.gz` file.

All that to say: `.igd` is a good format if you want fast access to variants in a relatively small file.

