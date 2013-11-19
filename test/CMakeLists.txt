########################################################################
# CMake build script for Google Test.
#
# To run the tests for Google Test itself on Linux, use 'make test' or
# ctest.  You can select which tests to run using 'ctest -R regex'.
# For more options, run 'ctest --help'.

# BUILD_SHARED_LIBS is a standard CMake variable, but we declare it here to
# make it prominent in the GUI.
option(BUILD_SHARED_LIBS "Build shared libraries (DLLs)." OFF)

# When other libraries are using a shared version of runtime libraries,
# Google Test also has to use one.
option(
  gtest_force_shared_crt
  "Use shared (DLL) run-time lib even when Google Test is built as static lib."
  OFF)

option(gtest_disable_pthreads "Disable uses of pthreads in gtest." OFF)

########################################################################
#
# Project-wide settings

# Name of the project.
#
# CMake files in this project can refer to the root source directory
# as ${gtest_SOURCE_DIR} and to the root binary directory as
# ${gtest_BINARY_DIR}.
# Language "C" is required for find_package(Threads).
#project(gtest CXX C)
#cmake_minimum_required(VERSION 2.6.2)

# Define helper functions and macros used by Google Test.
include(cmake/internal_utils.cmake)

config_compiler_and_linker()  # Defined in internal_utils.cmake.

#string(REGEX REPLACE "-ansi" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
string(REGEX REPLACE "-pedantic" " " CMAKE_CXX_FLAGS_DEVELOPER "${CMAKE_CXX_FLAGS_DEVELOPER}")
string(REGEX REPLACE "-pedantic" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")

# Where Google Test's .h files can be found.
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/gtest/gtest
  ${CMAKE_CURRENT_SOURCE_DIR}/gtest
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR})

# Where Google Test's libraries can be found.
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/gtest)

########################################################################
#
# Defines the gtest & gtest_main libraries.  User tests should link
# with one of them.

# Google Test libraries.  We build them using more strict warnings than what
# are used for other targets, to ensure that gtest can be compiled by a user
# aggressive about warnings.

#------------------------------------------------------------------------------
# Define test targets

cxx_library(gtest "${cxx_strict}" gtest/gtest-all.cc)
cxx_library(gtest_main "${cxx_strict}" gtest/gtest_main.cc)
target_link_libraries(gtest_main gtest)

set(GOSS_TEST_LIBS gtest_main ${PROJECT_NAME_LOWER} ${GOSS_TARGET_LINK_LIBRARIES})

#------------------------------------------------------------------------------
# Add all tests 

file(GLOB CPP_TESTS "*.cpp")

foreach (CPP_TEST_FULL_PATH ${CPP_TESTS})
  get_filename_component(CPP_TEST "${CPP_TEST_FULL_PATH}" NAME_WE)
  cxx_test(${CPP_TEST} "${GOSS_TEST_LIBS}")

endforeach()
