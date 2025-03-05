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
