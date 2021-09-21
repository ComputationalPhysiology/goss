cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(PYTHON_LIBRARY_DIR "/Users/henriknf/local/src/goss/venv/lib/python3.8/site-packages")
set(PYTHON_EXECUTABLE "/Users/henriknf/local/src/goss/venv/bin/python")


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(goss)

find_package(PythonLibs ${PY_VERSION} REQUIRED)
find_package(Boost 1.64.0 REQUIRED COMPONENTS system python)

include_directories("${CMAKE_SOURCE_DIR}/cpp/include/goss_bits")
include_directories("${CMAKE_SOURCE_DIR}/python")
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

file (GLOB SOURCE_FILES "cpp/src/*.cpp")
file (GLOB HEADER_FILES "cpp/include/goss_bits/*.h")
file (GLOB PYTHON_FILES "python/*.cpp" "python/*.hpp")

# Set up such that XCode organizes the files
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES} ${HEADER_FILES} ${PYTHON_FILES} )

include(pybind11.cmake)
pybind11_add_module(goss 
	${SOURCE_FILES}
	${HEADER_FILES}
	${PYTHON_FILES}
)

if(Boost_FOUND)
  target_link_libraries(goss PRIVATE ${Boost_LIBRARIES})
endif()



install(TARGETS goss
  COMPONENT python
  LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}"
  )