# CMake help file to add dependencies of goss

#------------------------------------------------------------------------------
# Run tests to find required packages

find_package(PythonInterp)
if (GOSS_ENABLE_PYTHON)
  find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT)

  # If Python is found, check for NumPy and SWIG
  if (PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND)
    find_package(NumPy REQUIRED)

    find_package(SWIG REQUIRED)
    if (${SWIG_VERSION} LESS 2.0)
      message(FATAL_ERROR " GOSS requires SWIG version 2.0 or greater. You have version ${SWIG_VERSION}. Set GOSS_ENABLE_PYTHON to False or install correct SWIG version.")
    endif()
    include(UseSWIG)
    set(PYTHON_FOUND TRUE)

    # Set numpy version define
    set(GOSS_PYTHON_DEFINITIONS -DNUMPY_VERSION_MAJOR=${NUMPY_VERSION_MAJOR} -DNUMPY_VERSION_MINOR=${NUMPY_VERSION_MINOR} -DNUMPY_VERSION_MICRO=${NUMPY_VERSION_MICRO})

    # Only set define for none depricated API for NUMPY version 1.7 and larger
    if(NUMPY_VERSION VERSION_GREATER 1.6.2)
      set(GOSS_PYTHON_DEFINITIONS ${DOLFIN_PYTHON_DEFINITIONS} -DNPY_NO_DEPRECATED_API=NPY_${NUMPY_VERSION_MAJOR}_${NUMPY_VERSION_MINOR}_API_VERSION)
    endif()

  endif()
endif()

# Look for Boost 
set(BOOST_ROOT $ENV{BOOST_DIR})
if (BOOST_ROOT)
  set(Boost_NO_SYSTEM_PATHS on)
endif()

# Prevent FindBoost.cmake from looking for system  Boost{foo}.cmake files
set(Boost_NO_BOOST_CMAKE true)

set(Boost_USE_MULTITHREADED $ENV{BOOST_USE_MULTITHREADED})
set(Boost_ADDITIONAL_VERSIONS 1.43 1.43.0 1.44 1.44.0 1.45 1.45.0 1.46 1.46.0 1.46.1 1.47 1.47.0 1.48 1.48.0 1.49 1.49.0 1.50 1.50.0)

find_package(Boost 1.36 QUIET REQUIRED)

set(GOSS_BOOST_COMPONENTS system thread)

find_package(Boost COMPONENTS ${GOSS_BOOST_COMPONENTS} REQUIRED)

# Check for OpenMP
if (GOSS_ENABLE_OPENMP)
  find_package(OpenMP)
  include(CheckOpenMP)
  check_openmp_unsigned_int_loop_control_variable(OPENMP_UINT_TEST_RUNS)
  if (NOT OPENMP_UINT_TEST_RUNS)
    set(OPENMP_FOUND FALSE)
  endif()
endif()

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
# Add include directories and libs of required packages
# Boost
list(APPEND GOSS_DEP_SYSTEM_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
list(APPEND GOSS_TARGET_LINK_LIBRARIES ${Boost_LIBRARIES})
list(APPEND APPEND GOSS_TARGET_LINK_LIBRARIES_DIRS ${Boost_LIBRARY_DIRS})

# OpenMP
if (GOSS_ENABLE_OPENMP AND OPENMP_FOUND)
  list(APPEND GOSS_CXX_DEFINITIONS "-DHAS_OPENMP")
  set(GOSS_CXX_FLAGS "${GOSS_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

  if (MINGW)
    list(APPEND GOSS_LINK_FLAGS "${OpenMP_CXX_FLAGS}")
  endif()

endif()

#------------------------------------------------------------------------------
# Set compiler flags, include directories and library dependencies

# Add compiler include directories
include_directories(${GOSS_SOURCE_DIR} ${GOSS_DEP_INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${GOSS_DEP_SYSTEM_INCLUDE_DIRECTORIES})

# Add CXX defintions
add_definitions(${GOSS_CXX_DEFINITIONS})
add_definitions(-DGOSS_VERSION="${GOSS_VERSION}")
add_definitions(${GOSS_PYTHON_DEFINITIONS})

# Add flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GOSS_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GOSS_LINK_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${GOSS_LINK_FLAGS}")

#------------------------------------------------------------------------------
# If everything needed for the Python extension module is found

if (PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND AND SWIG_FOUND AND Boost_FOUND)
    set(ENABLE_PYTHON_EXTENSION_MODULE TRUE)
endif()
