cmake_minimum_required(VERSION 3.16)

project(goss_pybind11 VERSION "0.4.2" LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(OpenMP)

if (TARGET goss)
  add_library(Goss::goss ALIAS goss)
else()
  # Find goss (C++)
  find_package(Goss REQUIRED)
endif()


if(NOT SKBUILD)
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.18.0")
    find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
  else()
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
  endif()
endif()

execute_process(
  COMMAND
    "${Python3_EXECUTABLE}" -c
    "import pybind11; print(pybind11.get_cmake_dir())"
  OUTPUT_VARIABLE _tmp_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")

find_package(pybind11 REQUIRED CONFIG HINTS ${PYBIND11_DIR} ${PYBIND11_ROOT}
  $ENV{PYBIND11_DIR} $ENV{PYBIND11_ROOT})


# Create the binding library
pybind11_add_module(_gosscpp MODULE wrapper.cpp)
target_link_libraries(_gosscpp PUBLIC Goss::goss)

# Add strict compiler flags
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-Wall -Werror  -Wextra -Wno-comment -pedantic" HAVE_PEDANTIC)

if(HAVE_PEDANTIC)
  target_compile_options(_gosscpp PRIVATE -Wall;-Wextra;-Werror;-Wno-comment;-pedantic)
endif()

# In Debug mode override pybind11 symbols visibility. Symbols must be
# visible to backtrace_symbols() to produce nice logs
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_options(_gosscpp PRIVATE "-fvisibility=default")
endif()

# scikit-build specific steps
if(SKBUILD)
  # Add relative rpath so _gosscpp can find the Goss::goss library
  # when the build is relocated
  if(APPLE)
    set_target_properties(_gosscpp PROPERTIES INSTALL_RPATH "@loader_path/lib")
  else()
    set_target_properties(_gosscpp PROPERTIES INSTALL_RPATH "$ORIGIN/lib")
  endif()

  install(TARGETS _gosscpp DESTINATION .)
endif()
#
