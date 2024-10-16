#!/bin/bash

set -ev

clang-format --dry-run -Werror igdtools/*.cpp
clang-format --dry-run -Werror examples/*.cpp
clang-format --dry-run -Werror *.hpp
