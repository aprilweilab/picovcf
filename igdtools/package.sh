#!/bin/bash

set -ev

cd /io/igdtools/

PKG=igdtools
VER=2

# Build the source distribution
/opt/python/cp38-cp38/bin/pip install setuptools
/opt/python/cp38-cp38/bin/python setup.py sdist

echo "Building for Python 3.8"
/opt/python/cp38-cp38/bin/python setup.py bdist_wheel
echo "Building for Python 3.9"
/opt/python/cp39-cp39/bin/pip install setuptools
/opt/python/cp39-cp39/bin/python setup.py bdist_wheel
echo "Building for Python 3.10"
/opt/python/cp310-cp310/bin/pip install setuptools
/opt/python/cp310-cp310/bin/python setup.py bdist_wheel
echo "Building for Python 3.11"
/opt/python/cp311-cp311/bin/pip install setuptools
/opt/python/cp311-cp311/bin/python setup.py bdist_wheel
echo "Building for Python 3.12"
/opt/python/cp312-cp312/bin/pip install setuptools
/opt/python/cp312-cp312/bin/python setup.py bdist_wheel
echo "Building for Python 3.13"
/opt/python/cp313-cp313/bin/pip install setuptools
/opt/python/cp313-cp313/bin/python setup.py bdist_wheel

cd /io/igdtools/dist
auditwheel repair --plat manylinux_2_24_x86_64 ${PKG}-${VER}.*-cp38-cp38-linux_x86_64.whl
auditwheel repair --plat manylinux_2_24_x86_64 ${PKG}-${VER}.*-cp39-cp39-linux_x86_64.whl
auditwheel repair --plat manylinux_2_24_x86_64 ${PKG}-${VER}.*-cp310-cp310-linux_x86_64.whl
auditwheel repair --plat manylinux_2_24_x86_64 ${PKG}-${VER}.*-cp311-cp311-linux_x86_64.whl
auditwheel repair --plat manylinux_2_24_x86_64 ${PKG}-${VER}.*-cp312-cp312-linux_x86_64.whl
auditwheel repair --plat manylinux_2_24_x86_64 ${PKG}-${VER}.*-cp313-cp313-linux_x86_64.whl
