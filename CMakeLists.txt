cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(cpp LANGUAGES CXX)

cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Use these to variables to set which python to use - should find a better way
# set(PYTHON_LIBRARY_DIR "/Users/henriknf/local/src/goss/venv/lib/python3.9/site-packages")
# set(PYTHON_EXECUTABLE "/Users/henriknf/local/src/goss/venv/bin/python")


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP REQUIRED)
find_package(PythonLibs ${PY_VERSION} REQUIRED)
# find_package(Boost 1.64.0 REQUIRED)
# find_package(Boost 1.64.0 REQUIRED COMPONENTS system python)
# find_package(Boost COMPONENTS python REQUIRED)


find_package(Boost 1.56 QUIET REQUIRED)

# Boost public/private libraries to link to.
# Note: These should all be private as they do not appear in the
# GOSS public interface , but there is a linking issues with older
# Boost or CMake. Ubuntu 16.04 requires linking GOSS programs with
# filesystem, whereas Ubuntu 16.10 and macOS (Homebrew) do not.
if (Boost_VERSION VERSION_LESS 106100)
  set(GOSS_BOOST_COMPONENTS_PUBLIC filesystem timer)
  set(GOSS_BOOST_COMPONENTS_PRIVATE program_options iostreams)
else()
  set(GOSS_BOOST_COMPONENTS_PUBLIC timer)
  set(GOSS_BOOST_COMPONENTS_PRIVATE filesystem program_options iostreams)
endif()

# Find required Boost libraries
find_package(Boost COMPONENTS ${GOSS_BOOST_COMPONENTS_PUBLIC} ${GOSS_BOOST_COMPONENTS_PRIVATE} REQUIRED)

# add_library(Boost::boost INTERFACE IMPORTED)
# set_target_properties(Boost::boost PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}")

include_directories("${CMAKE_SOURCE_DIR}/cpp/include/goss_bits")
include_directories("${CMAKE_SOURCE_DIR}/python")
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

file (GLOB SOURCES "cpp/src/*.cpp")
file (GLOB HEADERS "cpp/include/goss_bits/*.h")
file (GLOB BINDINGS "python/*.cpp" "python/*.hpp")

# Set up such that XCode organizes the files
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCES} ${HEADERS} ${BINDINGS} )

include(pybind11.cmake)
# pybind11_add_module(gos
# 	${SOURCES}
# 	${HEADERS}
# 	${BINDINGS}
# )

add_library(cpp MODULE "${SOURCES}" "${BINDINGS}")


target_link_libraries(cpp PRIVATE pybind11::module OpenMP::OpenMP_CXX Boost::boost ${PYTHON_LIBRARIES})
# target_include_directories(goss PRIVATE ${PYTHON_INCLUDE_DIRS})
set_target_properties(cpp PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                     SUFFIX "${PYTHON_MODULE_EXTENSION}")

# install(TARGETS goss
#   COMPONENT python
#   LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}"
#   )