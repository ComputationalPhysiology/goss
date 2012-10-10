# CMake help file to add dependencies of goss

#------------------------------------------------------------------------------
# Run tests to find required packages

find_package(PythonInterp)
if (GOSS_ENABLE_PYTHON)
  find_package(PythonLibs)

  # If Python is found, check for NumPy and SWIG
  if (PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND)
    find_package(NumPy REQUIRED)

    find_package(SWIG REQUIRED)
    if (${SWIG_VERSION} LESS 2.0)
      message(FATAL_ERROR " GOSS requires SWIG version 2.0 or greater. You have version ${SWIG_VERSION}. Set GOSS_ENABLE_PYTHON to False or install correct SWIG version.")
    endif()
    include(UseSWIG)
    set(PYTHON_FOUND TRUE)
  endif()
endif()

# Look for Boost (scoped_array support)
set(Boost_ADDITIONAL_VERSIONS 1.43 1.43.0 1.44 1.44.0 1.45 1.45.0 1.46 1.46.0 
    1.46.1 1.47 1.47.0 1.48 1.48.0 1.49 1.49.0 1.50 1.50.0)
set(BOOST_ROOT $ENV{BOOST_DIR})
find_package(Boost 1.36 QUIET)

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
# Add include directories and libs of required packages
# Boost
list(APPEND GOSS_DEP_SYSTEM_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})

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

# Add flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GOSS_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GOSS_LINK_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${GOSS_LINK_FLAGS}")

