cmake_minimum_required(VERSION 3.18.2)

project(goss_pybind11 VERSION "0.1.0"  LANGUAGES CXX)

include(GNUInstallDirs)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# See https://gitlab.kitware.com/cmake/cmake/-/issues/16414
if (TARGET goss)
  add_library(Goss::goss ALIAS goss)
else()
  # Find Goss (C++)
  find_package(Goss REQUIRED)
endif()

# find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
# find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

if(SKBUILD)
  # Scikit-Build does not add your site-packages to the search path
  # automatically, so we need to add it or the pybind11 specific
  # directory here.
  execute_process(
    COMMAND
      "${Python3_EXECUTABLE}" -c
      "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_VARIABLE _tmp_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
  list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
endif()

find_package(pybind11 REQUIRED CONFIG HINTS ${PYBIND11_DIR} ${PYBIND11_ROOT}
  $ENV{PYBIND11_DIR} $ENV{PYBIND11_ROOT})

set(python_module_name _gosscpp)

# Create the binding library
pybind11_add_module(_gosscpp MODULE wrapper.cpp)
target_link_libraries(_gosscpp PRIVATE Goss::goss)

# Add strict compiler flags
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-Wall -Werror  -Wextra -Wno-comment -pedantic" HAVE_PEDANTIC)
if (HAVE_PEDANTIC)
  target_compile_options(_gosscpp PRIVATE -Wall;-Wextra;-Werror;-Wno-comment;-pedantic)
endif()

# In Debug mode override pybind11 symbols visibility. Symbols must be
# visible to backtrace_symbols() to produce nice logs
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_options(_gosscpp PRIVATE "-fvisibility=default")
endif()

# scikit-build specific steps
if (SKBUILD)
  # Add relative rpath so _gosscpp can find the Goss::goss libray
  # when the build is relocated
  if (APPLE)
    set_target_properties(_gosscpp PROPERTIES INSTALL_RPATH "@loader_path/lib")
  else()
    set_target_properties(_gosscpp PROPERTIES INSTALL_RPATH "$ORIGIN/lib")
  endif()

  install(TARGETS _gosscpp DESTINATION .)
endif()
