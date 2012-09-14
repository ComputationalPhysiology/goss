# Top level CMakeLists.txt file for GOSS

# Require CMake 2.6
cmake_minimum_required(VERSION 2.6)

#------------------------------------------------------------------------------
# Set project name and version number

project(goss)

set(GOSS_VERSION_MAJOR "0")
set(GOSS_VERSION_MINOR "9")
set(GOSS_VERSION_MICRO "0")
set(GOSS_VERSION "${GOSS_VERSION_MAJOR}.${GOSS_VERSION_MINOR}.${GOSS_VERSION_MICRO}")

# Set special link option, see `cmake --help-policy CMP0003`
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif()

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)

# Set location of our FindFoo.cmake modules
set(GOSS_CMAKE_DIR "${GOSS_SOURCE_DIR}/cmake" CACHE INTERNAL "")
#set(CMAKE_MODULE_PATH "${GOSS_CMAKE_DIR}/modules")

#------------------------------------------------------------------------------
# Options

option(GOSS_ENABLE_PYTHON "Enable Python extensions." OFF)
option(GOSS_WITH_LIBRARY_VERSION "Build with library version information." ON)
option(GOSS_ENABLE_UNIT_TESTS "Build unit tests." ON)
option(BUILD_SHARED_LIBS "Build GOSS with shared libraries." ON)
option(CMAKE_SKIP_RPATH "Do not add runtime paths when using shared libraries." OFF)
option(CMAKE_INSTALL_RPATH_USE_LINK_PATH "Add paths to linker search and installed rpath." ON)
option(CMAKE_USE_RELATIVE_PATHS "Use relative paths in makefiles and projects." OFF)

#------------------------------------------------------------------------------
# Run tests to find required packages

find_package(PythonInterp QUIET)
find_package(PythonLibs QUIET)

# Find SWIG
find_package(SWIG QUIET)
include(UseSWIG)

# Look for Boost (shared_ptr support for Python wrappers)
#set(Boost_ADDITIONAL_VERSIONS 1.43 1.43.0)
#set(BOOST_ROOT $ENV{BOOST_DIR})
#find_package(Boost 1.36 QUIET)

#------------------------------------------------------------------------------
# Get installation path for Python modules

# Get Python module path from distutils
if (PYTHONINTERP_FOUND)

  if (NOT DEFINED GOSS_INSTALL_PYTHON_EXT_DIR)
    # Get path for platform-dependent Python modules (since we install a binary libary)
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))"
      OUTPUT_VARIABLE GOSS_INSTALL_PYTHON_EXT_DIR
      )
    # Strip off CMAKE_INSTALL_PREFIX (is added later by CMake)
    string(REGEX REPLACE "${CMAKE_INSTALL_PREFIX}(/|\\\\)([^ ]*)" "\\2"
      GOSS_INSTALL_PYTHON_EXT_DIR "${GOSS_INSTALL_PYTHON_EXT_DIR}")
    set(GOSS_INSTALL_PYTHON_EXT_DIR ${GOSS_INSTALL_PYTHON_EXT_DIR}
      CACHE PATH "Python extension module installation directory.")
  endif()

  if (NOT DEFINED GOSS_INSTALL_PYTHON_MODULE_DIR)
    # Get path for pure Python modules
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=False, prefix='${CMAKE_INSTALL_PREFIX}'))"
      OUTPUT_VARIABLE GOSS_INSTALL_PYTHON_MODULE_DIR
      )
    # Strip off CMAKE_INSTALL_PREFIX (is added later by CMake)
    string(REGEX REPLACE "${CMAKE_INSTALL_PREFIX}(/|\\\\)([^ ]*)" "\\2"
      GOSS_INSTALL_PYTHON_MODULE_DIR "${GOSS_INSTALL_PYTHON_MODULE_DIR}")
    set(GOSS_INSTALL_PYTHON_MODULE_DIR ${GOSS_INSTALL_PYTHON_MODULE_DIR}
      CACHE PATH "Python module installation directory.")
  endif()
endif (PYTHONINTERP_FOUND)

#------------------------------------------------------------------------------
# Build SWIG extension and install

if (GOSS_ENABLE_PYTHON AND SWIG_FOUND)# AND Boost_FOUND)

  # Check SWIG version
  if (${SWIG_VERSION} LESS 2.0)
      message(ERROR " GOSS requires SWIG version 2.0 or greater. You have version ${SWIG_VERSION}. Set GOSS_ENABLE_PYTHON to False or install correct SWIG version.")
    endif()

  # Default to release build (can be overridden by user)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
  endif()

  # Make build directory for SWIG-generated C++ file
  FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/src/ufc")

  # In CMake 2.6 PYTHON_INCLUDE_DIRS was named PYTHON_INCLUDE_PATH
  if (NOT DEFINED PYTHON_INCLUDE_DIRS)
    set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_PATH})
  endif()

  # Include directories
  include_directories(${GOSS_SOURCE_DIR}/goss/swig ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})

  # Set module name
  set(SWIG_MODULE_NAME goss)

  # SWIG flags
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
    )

  # SWIG sources
  set(SWIG_SOURCES goss/swig/goss.i)
  set(SWIG_MODULE_${SWIG_MODULE_NAME}_EXTRA_DEPS ${GOSS_H})
  set_source_files_properties(${SWIG_SOURCES} PROPERTIES CPLUSPLUS ON)
  swig_add_module(${SWIG_MODULE_NAME} python ${SWIG_SOURCES})

  # Is this required?
  swig_link_libraries(goss ${PYTHON_LIBRARIES})
  get_target_property(SWIG_MODULE_LOCATION ${SWIG_MODULE_goss_REAL_NAME} LOCATION)
  message("SWIG_MODULE_LOCATION ${SWIG_MODULE_LOCATION}")

  # Install the swig file
  install(FILES SWIG_SOURCES
    DESTINATION ${GOSS_INCLUDE_DIR}/swig
    COMPONENT Development
    )

  # Install _goss.so and goss.py
  install(FILES
    ${SWIG_MODULE_LOCATION}
    ${CMAKE_CURRENT_BINARY_DIR}/goss.py src/goss/__init__.py
    DESTINATION ${GOSS_INSTALL_PYTHON_EXT_DIR}/goss
    COMPONENT Development
    )

endif()

#------------------------------------------------------------------------------
# Installation of DOLFIN library

# Set GOSS install sub-directories
#set(GOSS_BIN_DIR "bin" CACHE PATH "Binary installation directory.")
set(GOSS_LIB_DIR "lib" CACHE PATH "Library installation directory.")
set(GOSS_INCLUDE_DIR "include" CACHE PATH "C/C++ header installation directory.")
set(GOSS_PKGCONFIG_DIR "lib/pkgconfig" CACHE PATH "pkg-config file installation directory.")
#set(GOSS_SHARE_DIR "share/goss" CACHE PATH "Shared data installation directory.")
#set(GOSS_MAN_DIR "share/man" CACHE PATH "Manual page installation directory.")
#set(GOSS_DOC_DIR "${GOSS_SHARE_DIR}/doc" CACHE PATH "GOSS Documentation directory.")

# Append the library version information to the library target properties
if (GOSS_WITH_LIBRARY_VERSION)
  string(REPLACE "+" "" GOSS_LIBRARY_VERSION ${GOSS_VERSION})
  # This setting of SOVERSION assumes that any API change
  # will increment either the minor or major version number.
  set(GOSS_LIBRARY_PROPERTIES ${GOSS_LIBRARY_PROPERTIES}
    VERSION ${GOSS_LIBRARY_VERSION}
    SOVERSION ${GOSS_VERSION_MAJOR}.${GOSS_VERSION_MINOR}
  )
endif()

# Add source directory
add_subdirectory(goss)

#------------------------------------------------------------------------------
# If enabling unit testing add test subdirectory

if (GOSS_ENABLE_UNIT_TESTS)
  add_subdirectory(test)
  
  # Add target "run_quicktest" for running only Python unit tests
  add_custom_target(run_tests
    COMMAND ode_test
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test")
endif()

#------------------------------------------------------------------------------
# Generate GOSSConfig.cmake file

#configure_file(${GOSS_CMAKE_DIR}/GOSSConfig.cmake.in GOSSConfig.cmake @ONLY)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GOSSConfig.cmake
#  DESTINATION ${GOSS_CMAKE_CONFIG_DIR}
#  COMPONENT Development
#  )

#------------------------------------------------------------------------------
# Generate GOSSConfigVersion.cmake file

#configure_file(${GOSS_CMAKE_DIR}/GOSSConfigVersion.cmake.in
#  GOSSConfigVersion.cmake @ONLY)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GOSSConfigVersion.cmake
#  DESTINATION ${GOSS_CMAKE_CONFIG_DIR}
#  COMPONENT Development
#  )

#------------------------------------------------------------------------------
# Generate pkg-config config file (goss.pc)

#configure_file(${GOSS_CMAKE_DIR}/goss.pc.in goss.pc @ONLY)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/goss.pc
#  DESTINATION ${GOSS_PKGCONFIG_DIR}
#  COMPONENT Development
#  )

#------------------------------------------------------------------------------

