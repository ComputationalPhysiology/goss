project(swig)

set(SWIG_MODULE_NAME cpp)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Remove '-Werror', '-Wall' and 'pedantic' flags (if present) when compiling SWIG-generated files
string(REGEX REPLACE "-Wall" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
string(REGEX REPLACE "-Wall" " " CMAKE_CXX_FLAGS_DEVELOPER "${CMAKE_CXX_FLAGS_DEVELOPER}")
string(REGEX REPLACE "-Werror=format-security" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
string(REGEX REPLACE "-Werror=format-security" " " CMAKE_CXX_FLAGS_DEVELOPER "${CMAKE_CXX_FLAGS_DEVELOPER}")
string(REGEX REPLACE "-Werror" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
string(REGEX REPLACE "-Werror" " " CMAKE_CXX_FLAGS_DEVELOPER "${CMAKE_CXX_FLAGS_DEVELOPER}")
string(REGEX REPLACE "-pedantic" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
string(REGEX REPLACE "-pedantic" " " CMAKE_CXX_FLAGS_DEVELOPER "${CMAKE_CXX_FLAGS_DEVELOPER}")

# Set common include directoris
include_directories(${GOSS_SOURCE_DIR} ${GOSS_SOURCE_DIR}/goss/swig ${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIR})

# Set SWIG flags 
set(CMAKE_SWIG_FLAGS
  -module ${SWIG_MODULE_NAME}
  -shadow
  -modern
  -modernargs
  -fastdispatch
  -fvirtual
  -nosafecstrings
  -noproxydel
  -fastproxy
  -fastinit
  -fastunpack
  -fastquery
  -nobuildnone
  -Iinclude/swig
  ${GOSS_CXX_DEFINITIONS}
  )

# Get all SWIG interface files

# All header files are dependencies
file(GLOB GOSS_HEADERS RELATIVE ${CMAKE_SOURCE_DIR} ../*.h)
file(GLOB GLOBAL_INTERFACE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.i)

list(APPEND GOSS_INTERFACE ${GOSS_HEADERS})
list(APPEND GOSS_INTERFACE ${GLOBAL_INTERFACE_FILES})

# This prevents swig being run unnecessarily
set_source_files_properties(goss.i PROPERTIES SWIG_MODULE_NAME ${SWIG_MODULE_NAME})

# Tell CMake SWIG has generated a C++ file
set_source_files_properties(goss.i PROPERTIES CPLUSPLUS ON)

# Generate SWIG files in 
set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR})

# Tell CMake which SWIG interface files should be checked for changes when recompile
set(SWIG_MODULE_${SWIG_MODULE_NAME}_EXTRA_DEPS ${GOSS_SWIG_DEPENDENCIES})

# Tell CMake to run SWIG on goss.i and to link against libgoss
swig_add_module(${SWIG_MODULE_NAME} python goss.i)
swig_link_libraries(${SWIG_MODULE_NAME} goss ${PYTHON_LIBRARIES})

# Install Python .py files
get_target_property(SWIG_MODULE_LOCATION ${SWIG_MODULE_${SWIG_MODULE_NAME}_REAL_NAME} LOCATION)

install(FILES
  ${SWIG_MODULE_LOCATION} ${CMAKE_CURRENT_BINARY_DIR}/${SWIG_MODULE_NAME}.py
  DESTINATION ${GOSS_INSTALL_PYTHON_MODULE_DIR}/goss
  COMPONENT RuntimeLibraries
  )

# Install all SWIG interface files
install(FILES ${GLOBAL_INTERFACE_FILES}
  DESTINATION 
  ${GOSS_INCLUDE_DIR}/goss/swig
  COMPONENT Development)

