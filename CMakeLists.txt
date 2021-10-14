cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(cpp LANGUAGES CXX)

# https://github.com/johnaparker/pybind_examples

cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(GOSS_LIB_DIR "${CMAKE_INSTALL_LIBDIR}" CACHE PATH "Library installation directory.")
set(GOSS_INCLUDE_DIR "${CMAKE_INSTALL_INCLUDEDIR}" CACHE PATH "C/C++ header installation directory.")
set(GOSS_PKGCONFIG_DIR "${CMAKE_INSTALL_LIBDIR}/pkgconfig" CACHE PATH "pkg-config file installation directory.")

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
# include_directories(${Boost_INCLUDE_DIRS})

file (GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/goss/*.cpp")
file (GLOB HEADERS "${PROJECT_SOURCE_DIR}/include/goss/*.h")
file (GLOB BINDINGS "${PROJECT_SOURCE_DIR}/src/python/goss/*.cpp" "${PROJECT_SOURCE_DIR}/src/python/goss/*.hpp")

include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/src
)
if(BUILD_TESTS)
    add_subdirectory("${CMAKE_SOURCE_DIR}/tests/cpp")
endif()

add_library(goss SHARED ${SOURCES} ${HEADERS})
target_link_libraries(goss PUBLIC OpenMP::OpenMP_CXX Boost::program_options Boost::system Boost::thread)
install(TARGETS goss
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# Build demo
# add_executable(oscilator oscilator.cpp)
# target_link_libraries(oscilator goss)
add_executable(main main.cpp)
target_link_libraries(main goss)

# # Build python module
include(pybind11.cmake)
pybind11_add_module(cpp ${BINDINGS})
target_link_libraries(cpp PRIVATE pybind11::module goss OpenMP::OpenMP_CXX Boost::program_options Boost::system Boost::thread ${PYTHON_LIBRARIES})
set_target_properties(cpp PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                     SUFFIX "${PYTHON_MODULE_EXTENSION}")