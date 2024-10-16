#!/bin/bash

set -ev

clang-format -i igdtools/*.cpp
clang-format -i examples/*.cpp
clang-format -i *.hpp
