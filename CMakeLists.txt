# Top level CMakeLists.txt file for GOSSOL

# Require CMake 2.6
cmake_minimum_required(VERSION 2.6)

#------------------------------------------------------------------------------
# Set project name and version number

project(GOSSOL)
set(GOSSOL_VERSION_MAJOR "0")
set(GOSSOL_VERSION_MINOR "9")
set(GOSSOL_VERSION_MICRO "0")
set(GOSSOL_VERSION "${GOSSOL_VERSION_MAJOR}.${GOSSOL_VERSION_MINOR}.${GOSSOL_VERSION_MICRO}")

# Set special link option, see `cmake --help-policy CMP0003`
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif()

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)

# Set location of our FindFoo.cmake modules
set(GOSSOL_CMAKE_DIR "${GOSSOL_SOURCE_DIR}/cmake" CACHE INTERNAL "")
#set(CMAKE_MODULE_PATH "${GOSSOL_CMAKE_DIR}/modules")

#------------------------------------------------------------------------------
# Options

option(GOSSOL_ENABLE_PYTHON "Enable Python extensions." OFF)
option(GOSSOL_WITH_LIBRARY_VERSION "Build with library version information." ON)
option(BUILD_SHARED_LIBS "Build GOSSOL with shared libraries." ON)
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

  if (NOT DEFINED GOSSOL_INSTALL_PYTHON_EXT_DIR)
    # Get path for platform-dependent Python modules (since we install a binary libary)
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))"
      OUTPUT_VARIABLE GOSSOL_INSTALL_PYTHON_EXT_DIR
      )
    # Strip off CMAKE_INSTALL_PREFIX (is added later by CMake)
    string(REGEX REPLACE "${CMAKE_INSTALL_PREFIX}(/|\\\\)([^ ]*)" "\\2"
      GOSSOL_INSTALL_PYTHON_EXT_DIR "${GOSSOL_INSTALL_PYTHON_EXT_DIR}")
    set(GOSSOL_INSTALL_PYTHON_EXT_DIR ${GOSSOL_INSTALL_PYTHON_EXT_DIR}
      CACHE PATH "Python extension module installation directory.")
  endif()

  if (NOT DEFINED GOSSOL_INSTALL_PYTHON_MODULE_DIR)
    # Get path for pure Python modules
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=False, prefix='${CMAKE_INSTALL_PREFIX}'))"
      OUTPUT_VARIABLE GOSSOL_INSTALL_PYTHON_MODULE_DIR
      )
    # Strip off CMAKE_INSTALL_PREFIX (is added later by CMake)
    string(REGEX REPLACE "${CMAKE_INSTALL_PREFIX}(/|\\\\)([^ ]*)" "\\2"
      GOSSOL_INSTALL_PYTHON_MODULE_DIR "${GOSSOL_INSTALL_PYTHON_MODULE_DIR}")
    set(GOSSOL_INSTALL_PYTHON_MODULE_DIR ${GOSSOL_INSTALL_PYTHON_MODULE_DIR}
      CACHE PATH "Python module installation directory.")
  endif()
endif (PYTHONINTERP_FOUND)

#------------------------------------------------------------------------------
# Build SWIG extension and install

if (GOSSOL_ENABLE_PYTHON AND SWIG_FOUND)# AND Boost_FOUND)

  # Check SWIG version
  if (${SWIG_VERSION} LESS 2.0)
      message(ERROR " GOSSOL requires SWIG version 2.0 or greater. You have version ${SWIG_VERSION}. Set GOSSOL_ENABLE_PYTHON to False or install correct SWIG version.")
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
  include_directories(${GOSSOL_SOURCE_DIR}/gossol/swig ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})

  # Set module name
  set(SWIG_MODULE_NAME gossol)

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
  set(SWIG_SOURCES gossol/swig/gossol.i)
  set(SWIG_MODULE_${SWIG_MODULE_NAME}_EXTRA_DEPS ${GOSSOL_H})
  set_source_files_properties(${SWIG_SOURCES} PROPERTIES CPLUSPLUS ON)
  swig_add_module(${SWIG_MODULE_NAME} python ${SWIG_SOURCES})

  # Is this required?
  swig_link_libraries(gossol ${PYTHON_LIBRARIES})
  get_target_property(SWIG_MODULE_LOCATION ${SWIG_MODULE_gossol_REAL_NAME} LOCATION)
  message("SWIG_MODULE_LOCATION ${SWIG_MODULE_LOCATION}")

  # Install the swig file
  install(FILES SWIG_SOURCES
    DESTINATION ${GOSSOL_INCLUDE_DIR}/swig
    COMPONENT Development
    )

  # Install _gossol.so and gossol.py
  install(FILES
    ${SWIG_MODULE_LOCATION}
    ${CMAKE_CURRENT_BINARY_DIR}/gossol.py src/gossol/__init__.py
    DESTINATION ${GOSSOL_INSTALL_PYTHON_EXT_DIR}/gossol
    COMPONENT Development
    )

endif()

#------------------------------------------------------------------------------
# Installation of DOLFIN library

# Set GOSSOL install sub-directories
#set(GOSSOL_BIN_DIR "bin" CACHE PATH "Binary installation directory.")
set(GOSSOL_LIB_DIR "lib" CACHE PATH "Library installation directory.")
set(GOSSOL_INCLUDE_DIR "include" CACHE PATH "C/C++ header installation directory.")
set(GOSSOL_PKGCONFIG_DIR "lib/pkgconfig" CACHE PATH "pkg-config file installation directory.")
#set(GOSSOL_SHARE_DIR "share/gossol" CACHE PATH "Shared data installation directory.")
#set(GOSSOL_MAN_DIR "share/man" CACHE PATH "Manual page installation directory.")
#set(GOSSOL_DOC_DIR "${GOSSOL_SHARE_DIR}/doc" CACHE PATH "GOSSOL Documentation directory.")

# Append the library version information to the library target properties
if (GOSSOL_WITH_LIBRARY_VERSION)
  string(REPLACE "+" "" GOSSOL_LIBRARY_VERSION ${GOSSOL_VERSION})
  # This setting of SOVERSION assumes that any API change
  # will increment either the minor or major version number.
  set(GOSSOL_LIBRARY_PROPERTIES ${GOSSOL_LIBRARY_PROPERTIES}
    VERSION ${GOSSOL_LIBRARY_VERSION}
    SOVERSION ${GOSSOL_VERSION_MAJOR}.${GOSSOL_VERSION_MINOR}
  )
endif()

# Add source directory
add_subdirectory(gossol)

#------------------------------------------------------------------------------
# Generate GOSSOLConfig.cmake file

#configure_file(${GOSSOL_CMAKE_DIR}/GOSSOLConfig.cmake.in GOSSOLConfig.cmake @ONLY)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GOSSOLConfig.cmake
#  DESTINATION ${GOSSOL_CMAKE_CONFIG_DIR}
#  COMPONENT Development
#  )

#------------------------------------------------------------------------------
# Generate GOSSOLConfigVersion.cmake file

#configure_file(${GOSSOL_CMAKE_DIR}/GOSSOLConfigVersion.cmake.in
#  GOSSOLConfigVersion.cmake @ONLY)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GOSSOLConfigVersion.cmake
#  DESTINATION ${GOSSOL_CMAKE_CONFIG_DIR}
#  COMPONENT Development
#  )

#------------------------------------------------------------------------------
# Generate pkg-config config file (gossol.pc)

#configure_file(${GOSSOL_CMAKE_DIR}/gossol.pc.in gossol.pc @ONLY)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gossol.pc
#  DESTINATION ${GOSSOL_PKGCONFIG_DIR}
#  COMPONENT Development
#  )

#------------------------------------------------------------------------------

