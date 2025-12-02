![](https://github.com/aprilweilab/picovcf/actions/workflows/cmake-multi-platform.yml/badge.svg)
![](https://readthedocs.org/projects/picovcf/badge/?version=latest)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/igdtools/README.html)

# picovcf

Single-header C++ library for fast/low-memory VCF (Variant Call Format) parsing. Gzipped VCF (.vcf.gz) is optionally supported.

**NEW!** Tabix indexes are now supported. See the `lowerBoundPosition()` method in the [documentation](https://picovcf.readthedocs.io/) for an example of how Tabix indexes are used by `picovcf`.

There are a lot of great tools for processing VCF files out there, but not many C++ libraries that are small (only parsing, no extra functionality) and easy to use. `picovcf` attempts to fill this niche by providing a header-only library using modern C++ (C++11) that allows clients to be selective about which parts of the VCF file get parsed.

Features:
* Fast and easy to use VCF(.GZ) parsing.
* Convert VCF(.GZ) to Indexable Genotype Data (IGD) format, which is a very simple format that is **more than 3x smaller than VCF.GZ at Biobank scale** and **more than 15x faster to read**
* Fast and easy to use IGD parsing.

More details on IGD can be found in our [preprint IGD paper](https://www.biorxiv.org/content/10.1101/2025.02.05.636549v1.abstract).

See also [pyigd](https://github.com/aprilweilab/pyigd/) if you want Python access to IGD files.

## Using the library

Either copy the latest header file (`picovcf.hpp`) into your project directly, or make use of something like git submodules to include https://github.com/aprilweilab/picovcf.

See the [vcfpp.cpp](https://github.com/aprilweilab/picovcf/blob/main/examples/vcfpp.cpp) for an example of how to use the APIs. Read [the docs](https://picovcf.readthedocs.io/en/latest/) for an overview of the API.

When building code that uses `picovcf.hpp`, define `VCF_GZ_SUPPORT=1` (`-DVCF_GZ_SUPPORT=1` on most compiler command lines) to enable zlib support for compressed VCF files.

## Build and run the tests/tools

picovcf does not need to be built to be used, since it is a single header that gets built as part of your project. However, if you want to build the tests and tools:

```
cd picovcf
mkdir build && cd build
cmake .. -DENABLE_VCF_GZ=ON
make
```

**NOTE**: `-DENABLE_VCF_GZ=ON` is optional, and links against `libz` in case you want to support `.vcf.gz` (compressed) files in the tools.

To convert from a `.vcf` or `.vcf.gz` file to `.igd`, run:
```
./igdtools <vcf filename> -o <output IGD filename>
```

Run `./igdtools --help` to see the full list of options. Here are some common tasks you might want to perform, besides VCF conversion:
* Pipe allele frequencies to a file: `./igdtools <input IGD> -a > allele.freq.tsv`
* View variant/sample statistics and header info: `./igdtools <input IGD> --stats --info`
* To, e.g., restrict to variants in base-pair range 10000,20000 add argument `--range 10000-20000`
* To restrict to variants with frequencies >=0.01: `--frange 0.01-1.0`
* Copy from one IGD to another: `./igdtools <input IGD> -o <output IGD>`
  * Only include variants in a certain range and with frequency: `./igdtools <input IGD> -o <output IGD> --range 100000-500000 --frange 0.01-0.5`

Finally, to run the unit tests:
```
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
* Simple format: all variants are expanded into binary variants. So if a Variant has `N` alternate alleles, then IGD will store that as `N` rows containing `0` (reference allele) or `1` (alternate allele). Each of these binary variants is stored as either a bitvector (non-sparse) or a list of sample indexes (sparse). A flag in the index indicates which way each variant is stored.
* Very small. Oftentimes smaller than compressed formats like `.vcf.gz` or `.bgen`. The more low-frequency mutations (such as for really large sample sizes) the smaller the file, assuming you are using the default implementation of dynamically choosing between sparse/non-sparse representation.

For example, the following are from chromosome 22 of a real dataset:
* `.vcf`: 11GB
* `.vcf.gz`: 203MB
* `.bgen`: 256MB
* `.igd`: 183MB

Converting the `.vcf.gz` to `.bgen` (via qctool) took 23 minutes, but converting to `.igd` only took 3 minutes. Furthermore, iteratively accessing all the variants (and genotype data) in the `.igd` file was approximately `15x` faster than accessing the same data in the `.vcf.gz` file (using `picovcf`). On Biobank-scale real datasets, IGD is on average 3.5x smaller than `.vcf.gz`.

### Installing just `igdtools`

You can install `igdtools` via PyPi: `pip install igdtools`

Or via [the conda package](https://bioconda.github.io/recipes/igdtools/README.html), after adding the [bioconda](https://bioconda.github.io/) channel: `conda install igdtools`

### How do I use IGD in my project?

* Clone [picovcf](https://github.com/aprilweilab/picovcf) and follow the instructions in this README to build the example tools for that library.
  * If you want to be able to convert `.vcf.gz` (compressed VCF) to IGD, make sure you build with `-DENABLE_VCF_GZ=ON`
* Use `igdtools` to convert and process files
* Do one of the following:
  * If your project is C++, copy [picovcf.hpp](https://github.com/aprilweilab/picovcf/blob/main/picovcf.hpp) into your project, `#include` it somewhere and then use according to the [documentation](https://picovcf.readthedocs.io/en/latest/)
  * If your project is Python, you can install [pyigd](https://github.com/aprilweilab/pyigd/) via `pip install pyigd` (see [the docs](https://pyigd.readthedocs.io/en/latest/))
