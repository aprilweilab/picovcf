# This file does not build picovcf (which is a single-header library, and does not need to
# be built separately from your project). It just builds tests.
cmake_minimum_required(VERSION 3.10)
project(picovcf)

option(ENABLE_VCF_GZ "Enable support for .vcf.gz (via ZLIB)" OFF)
option(FUZZING_SUPPORT "Build for fuzzing; only abort on true failures, not invalid input data" OFF)
option(BUILD_EXAMPLES "Build example tools" ON)

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

set(PICOVCF_TEST_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/test/test_helpers.cpp
    ${CMAKE_CURRENT_LIST_DIR}/test/test_igd.cpp
    ${CMAKE_CURRENT_LIST_DIR}/test/test_main.cpp
    ${CMAKE_CURRENT_LIST_DIR}/test/test_rotate.cpp
    ${CMAKE_CURRENT_LIST_DIR}/test/test_valid_vcf.cpp
  )

if (${ENABLE_VCF_GZ})
  add_compile_options(-DVCF_GZ_SUPPORT=1)
  set(LIBRARIES_TO_LINK z)
else()
  set(LIBRARIES_TO_LINK)
endif()

if (${FUZZING_SUPPORT})
  add_compile_options(-DFUZZING=1)
endif()

include_directories(${CMAKE_CURRENT_LIST_DIR})

# https://github.com/google/googletest/tree/master/googletest#incorporating-into-an-existing-cmake-project

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(CTest)
add_executable(picovcf_test ${PICOVCF_TEST_SOURCES})
target_link_libraries(picovcf_test gtest_main ${LIBRARIES_TO_LINK})
add_test(NAME picovcf_test COMMAND picovcf_test)

include(GoogleTest)
gtest_discover_tests(picovcf_test)

if(${BUILD_EXAMPLES})
  # vcfpp utility
  add_executable(vcfpp examples/vcfpp.cpp)
  target_link_libraries(vcfpp ${LIBRARIES_TO_LINK})
  # vcfconv utility
  add_executable(vcfconv examples/vcfconv.cpp)
  target_link_libraries(vcfconv ${LIBRARIES_TO_LINK})
  # igdpp utility
  add_executable(igdpp examples/igdpp.cpp)
  target_link_libraries(igdpp ${LIBRARIES_TO_LINK})
  # igdcp utility
  add_executable(igdcp examples/igdcp.cpp)
  target_link_libraries(igdcp ${LIBRARIES_TO_LINK})
  if (${ENABLE_VCF_GZ})
    # gzcat utility
    add_executable(gzcat examples/gzcat.cpp)
    target_link_libraries(gzcat ${LIBRARIES_TO_LINK})
  endif()
endif()
