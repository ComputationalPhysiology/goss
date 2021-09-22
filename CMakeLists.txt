cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(cpp LANGUAGES CXX)

cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

option(BUILD_TESTS "Build the test suite." ON)

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
find_package(Boost REQUIRED)

# find_package(Boost 1.56 QUIET REQUIRED)
# set(GOSS_BOOST_COMPONENTS_PUBLIC timer)
# set(GOSS_BOOST_COMPONENTS_PRIVATE filesystem program_options iostreams)

# # Find required Boost libraries
# find_package(Boost COMPONENTS ${GOSS_BOOST_COMPONENTS_PUBLIC} ${GOSS_BOOST_COMPONENTS_PRIVATE} REQUIRED)

# add_library(Boost::boost INTERFACE IMPORTED)
# set_target_properties(Boost::boost PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}")

# include_directories("${CMAKE_SOURCE_DIR}/cpp/include/goss_bits")
# include_directories("${CMAKE_SOURCE_DIR}/python")
# include_directories(${PYTHON_INCLUDE_DIRS})
# include_directories(${Boost_INCLUDE_DIRS})

file (GLOB SOURCES "src/goss/*.cpp")
file (GLOB HEADERS "include/goss/*.h")
file (GLOB BINDINGS "python/*.cpp" "python/*.hpp")

# Set up such that XCode organizes the files
# source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCES} ${HEADERS} ${BINDINGS} )

include(pybind11.cmake)
# pybind11_add_module(goss
# 	${SOURCES}
# 	${HEADERS}
# 	${BINDINGS}
# )


include_directories(
		${PROJECT_SOURCE_DIR}/include
		${PROJECT_SOURCE_DIR}/src
)
if(BUILD_TESTS)
    add_subdirectory("${CMAKE_SOURCE_DIR}/cpp/tests")
endif()

add_library(goss "${SOURCES}" "${HEADERS}")
target_link_libraries(goss PRIVATE OpenMP::OpenMP_CXX Boost::boost)

add_library(gosscpp MODULE "${SOURCES}" "${BINDINGS}" "${HEADERS}")
target_link_libraries(gosscpp PRIVATE pybind11::module OpenMP::OpenMP_CXX Boost::boost ${PYTHON_LIBRARIES})
# target_include_directories(goss PRIVATE ${PYTHON_INCLUDE_DIRS})
set_target_properties(gosscpp PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                     SUFFIX "${PYTHON_MODULE_EXTENSION}")

                                     
# install(TARGETS goss
#   COMPONENT python
#   LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}"
#   )