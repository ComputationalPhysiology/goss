# Top level CMakeLists.txt file for GOSOL

# Require CMake 2.6
cmake_minimum_required(VERSION 2.6)

#------------------------------------------------------------------------------
# Set project name and version number

project(GOSOL)
set(GOSOL_VERSION_MAJOR "0")
set(GOSOL_VERSION_MINOR "9")
set(GOSOL_VERSION_MICRO "0")
set(GOSOL_VERSION_STRING "${GOSOL_VERSION_MAJOR}.${GOSOL_VERSION_MINOR}.${GOSOL_VERSION_MICRO}")

# Set special link option, see `cmake --help-policy CMP0003`
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif()

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)

# Set location of our FindFoo.cmake modules
set(GOSOL_CMAKE_DIR "${GOSOL_SOURCE_DIR}/cmake" CACHE INTERNAL "")
#set(CMAKE_MODULE_PATH "${GOSOL_CMAKE_DIR}/modules")

#------------------------------------------------------------------------------
# Options

option(GOSOL_ENABLE_PYTHON "Enable Python extensions." OFF)

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

  if (NOT DEFINED GOSOL_INSTALL_PYTHON_EXT_DIR)
    # Get path for platform-dependent Python modules (since we install a binary libary)
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))"
      OUTPUT_VARIABLE GOSOL_INSTALL_PYTHON_EXT_DIR
      )
    # Strip off CMAKE_INSTALL_PREFIX (is added later by CMake)
    string(REGEX REPLACE "${CMAKE_INSTALL_PREFIX}(/|\\\\)([^ ]*)" "\\2"
      GOSOL_INSTALL_PYTHON_EXT_DIR "${GOSOL_INSTALL_PYTHON_EXT_DIR}")
    set(GOSOL_INSTALL_PYTHON_EXT_DIR ${GOSOL_INSTALL_PYTHON_EXT_DIR}
      CACHE PATH "Python extension module installation directory.")
  endif()

  if (NOT DEFINED GOSOL_INSTALL_PYTHON_MODULE_DIR)
    # Get path for pure Python modules
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=False, prefix='${CMAKE_INSTALL_PREFIX}'))"
      OUTPUT_VARIABLE GOSOL_INSTALL_PYTHON_MODULE_DIR
      )
    # Strip off CMAKE_INSTALL_PREFIX (is added later by CMake)
    string(REGEX REPLACE "${CMAKE_INSTALL_PREFIX}(/|\\\\)([^ ]*)" "\\2"
      GOSOL_INSTALL_PYTHON_MODULE_DIR "${GOSOL_INSTALL_PYTHON_MODULE_DIR}")
    set(GOSOL_INSTALL_PYTHON_MODULE_DIR ${GOSOL_INSTALL_PYTHON_MODULE_DIR}
      CACHE PATH "Python module installation directory.")
  endif()
endif (PYTHONINTERP_FOUND)

#------------------------------------------------------------------------------
# Target names and installation directories

# Set GOSOL install sub-directories
set(GOSOL_LIB_DIR "lib")
set(GOSOL_INCLUDE_DIR "include")
set(GOSOL_PKGCONFIG_DIR "lib/pkgconfig")
set(GOSOL_SHARE_DIR "share")
set(GOSOL_MAN_DIR "share/man")
set(GOSOL_CMAKE_CONFIG_DIR "share/ufc")

#------------------------------------------------------------------------------
# Build SWIG extension and install

if (GOSOL_ENABLE_PYTHON AND SWIG_FOUND)# AND Boost_FOUND)

  # Check SWIG version
  if (${SWIG_VERSION} LESS 2.0)
      message(ERROR " GOSOL requires SWIG version 2.0 or greater. You have version ${SWIG_VERSION}. Set GOSOL_ENABLE_PYTHON to False or install correct SWIG version.")
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
  include_directories(${GOSOL_SOURCE_DIR}/gosol/swig ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})

  # Set module name
  set(SWIG_MODULE_NAME gosol)

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
  set(SWIG_SOURCES gosol/swig/gosol.i)
  set(SWIG_MODULE_${SWIG_MODULE_NAME}_EXTRA_DEPS ${GOSOL_H})
  set_source_files_properties(${SWIG_SOURCES} PROPERTIES CPLUSPLUS ON)
  swig_add_module(${SWIG_MODULE_NAME} python ${SWIG_SOURCES})

  # Is this required?
  swig_link_libraries(gosol ${PYTHON_LIBRARIES})
  get_target_property(SWIG_MODULE_LOCATION ${SWIG_MODULE_gosol_REAL_NAME} LOCATION)
  message("SWIG_MODULE_LOCATION ${SWIG_MODULE_LOCATION}")

  # Install the swig file
  install(FILES SWIG_SOURCES
    DESTINATION ${GOSOL_INCLUDE_DIR}/swig
    COMPONENT Development
    )

  # Install _gosol.so and gosol.py
  install(FILES
    ${SWIG_MODULE_LOCATION}
    ${CMAKE_CURRENT_BINARY_DIR}/gosol.py src/gosol/__init__.py
    DESTINATION ${GOSOL_INSTALL_PYTHON_EXT_DIR}/gosol
    COMPONENT Development
    )

endif()

#------------------------------------------------------------------------------
# Add source directory
add_subdirectory(gosol)

#------------------------------------------------------------------------------
# Generate GOSOLConfig.cmake file

#configure_file(${GOSOL_CMAKE_DIR}/GOSOLConfig.cmake.in GOSOLConfig.cmake @ONLY)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GOSOLConfig.cmake
#  DESTINATION ${GOSOL_CMAKE_CONFIG_DIR}
#  COMPONENT Development
#  )

#------------------------------------------------------------------------------
# Generate GOSOLConfigVersion.cmake file

#configure_file(${GOSOL_CMAKE_DIR}/GOSOLConfigVersion.cmake.in
#  GOSOLConfigVersion.cmake @ONLY)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GOSOLConfigVersion.cmake
#  DESTINATION ${GOSOL_CMAKE_CONFIG_DIR}
#  COMPONENT Development
#  )

#------------------------------------------------------------------------------
# Generate pkg-config config file (gosol.pc)

#configure_file(${GOSOL_CMAKE_DIR}/gosol.pc.in gosol.pc @ONLY)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gosol.pc
#  DESTINATION ${GOSOL_PKGCONFIG_DIR}
#  COMPONENT Development
#  )

#------------------------------------------------------------------------------

