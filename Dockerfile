# Example to build:
#   docker build . -t picovcf:latest
#
# This file builds the picovcf examples, tests, and documentation. It mostly just serves as an example
# and an easy way to try things out. When you incorporate picovcf in your project, you only need the
# picovcf.hpp header file.
FROM ubuntu:22.04

RUN apt update
# For the base build
RUN apt install -y git build-essential cmake
# For build the documentation
RUN apt install -y python3 python3-setuptools python3-pip doxygen
RUN pip3 install wheel sphinx breathe sphinx-rtd-theme
# For .vcf.gz support
RUN apt install -y zlib1g-dev

COPY . /picovcf_src
RUN cd /picovcf_src \
    && mkdir docker_build \
    && cd docker_build \
    && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_VCF_GZ=ON \
    && make -j \
    && DOC_BUILD_DIR=$PWD sphinx-build -c ../doc/ -b html -Dbreathe_projects.picovcf=$PWD/doc/xml ../doc/ $PWD/doc/sphinx/
