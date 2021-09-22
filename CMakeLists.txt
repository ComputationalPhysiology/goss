cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(cpp LANGUAGES CXX)

# https://github.com/johnaparker/pybind_examples

cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

option(BUILD_TESTS "Build the test suite." OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP REQUIRED)
find_package(PythonLibs ${PY_VERSION} REQUIRED)
find_package(Boost COMPONENTS system thread program_options REQUIRED)


include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

file (GLOB SOURCES "src/goss/*.cpp")
# file (GLOB HEADERS "include/goss/*.h")
file (GLOB BINDINGS "src/python/goss/*.cpp" "src/python/goss/*.hpp")

include_directories(${PROJECT_SOURCE_DIR}/include)
if(BUILD_TESTS)
    add_subdirectory("${CMAKE_SOURCE_DIR}/tests/cpp")
endif()

# Build library
add_library(goss ${SOURCES})
target_link_libraries(goss PRIVATE OpenMP::OpenMP_CXX ${Boost_LIBRARIES})

# Build demo
add_executable(main main.cpp ${SOURCES})
target_link_libraries(main OpenMP::OpenMP_CXX ${Boost_LIBRARIES})

# Build python module
include(pybind11.cmake)
add_library(gosscpp MODULE ${SOURCES} ${BINDINGS})
target_link_libraries(gosscpp PRIVATE pybind11::module OpenMP::OpenMP_CXX ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
set_target_properties(gosscpp PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                     SUFFIX "${PYTHON_MODULE_EXTENSION}")

                                     
# install(TARGETS goss
#   COMPONENT python
#   LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}"
#   )