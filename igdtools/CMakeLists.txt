# This file does not build picovcf (which is a single-header library, and does not need to
# be built separately from your project). It just builds tests.
cmake_minimum_required(VERSION 3.10)
project(picovcf)

option(ENABLE_VCF_GZ "Enable support for .vcf.gz (via ZLIB)" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Check for optional POSIX functions
include(CheckSymbolExists)
check_symbol_exists(posix_fadvise "fcntl.h" HAVE_POSIX_FADVISE)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "In-source builds not allowed. Please make a build directory and run CMake from there.\n")
endif()

if (${ENABLE_VCF_GZ})
  add_compile_options(-DVCF_GZ_SUPPORT=1)
  set(LIBRARIES_TO_LINK z)
else()
  set(LIBRARIES_TO_LINK)
endif()

include_directories(${CMAKE_CURRENT_LIST_DIR})

add_executable(igdtools igdtools.cpp)
target_link_libraries(igdtools ${LIBRARIES_TO_LINK})
