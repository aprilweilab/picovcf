# Releasing new igdtools versions

This outlines the steps for releasing a new version of igdtools. This is just for the wrapper
that posts igdtools to PyPi, there are other ways to install (build manually, or use the
Docker image).

## Version numbering

Version numbering follows picovcf's version, always. Bump the version in setup.py whenever
a release is done on GitHub.

## Packaging for PyPi

Build the package distributions for PyPi. We build a source dist and then Linux binary distributions for recent Python versions. The container is from the [manylinux](https://github.com/pypa/manylinux) project.

```
# Remove dist/ to start fresh
rm -rf dist

# Pull the container for manylinux:
docker pull quay.io/pypa/manylinux_2_28_x86_64

# Run the packaging inside the container
docker run -v $PWD/../:/io -it quay.io/pypa/manylinux_2_28_x86_64 /io/igdtools/package.sh

# Fix file permissions from Docker
sudo chown -R ddehaas dist/
sudo chgrp -R ddehaas dist/

# Copy the source wheel to wheelhouse
cp ./dist/*.tar.gz ./dist/wheelhouse/

```

Test the source distribution:
```
pip install --force-reinstall ./dist/wheelhouse/igdtools-*.tar.gz
```

To upload to PyPi:
```
python3 -m twine upload dist/wheelhouse/*
```
