cmake_minimum_required(VERSION 3.16)

# Set the version
project(Goss VERSION "0.1.0" LANGUAGES CXX)
include(GNUInstallDirs)

if(SKBUILD)
  # Always use lib/ in the Python root. Otherwise, environment used for
  # wheels might place in lib64/, which CMake on some systems will not
  # find.
  set(CMAKE_INSTALL_LIBDIR lib)
endif()

include(FeatureSummary)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(APPLE)
  set(CMAKE_MACOSX_RPATH ON)
endif(APPLE)

# Options
option(BUILD_SHARED_LIBS "Build Goss with shared libraries." ON)
add_feature_info(BUILD_SHARED_LIBS BUILD_SHARED_LIBS "Build Goss with shared libraries.")
option(BUILD_TESTS "Build the test suite." OFF)

# Check if OpenMP is turned Off, otherwise set it to true
if(DEFINED ENV{OPENMP})
  set(USE_OPENMP $ENV{OPENMP})
else()
  set(USE_OPENMP 1)
endif()

if(USE_OPENMP)
  message("Try to find OpenMP")

  if(APPLE)
    find_package(OpenMP)

    if(NOT OpenMP_FOUND)
      # libomp 15.0+ from brew is keg-only, so have to search in other locations.
      # See https://github.com/Homebrew/homebrew-core/issues/112107#issuecomment-1278042927.
      execute_process(COMMAND brew --prefix libomp
        OUTPUT_VARIABLE HOMEBREW_LIBOMP_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE)
      set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I${HOMEBREW_LIBOMP_PREFIX}/include")
      set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${HOMEBREW_LIBOMP_PREFIX}/include")
      set(OpenMP_C_LIB_NAMES omp)
      set(OpenMP_CXX_LIB_NAMES omp)
      set(OpenMP_omp_LIBRARY ${HOMEBREW_LIBOMP_PREFIX}/lib/libomp.dylib)
      find_package(OpenMP)
    endif()
  else()
    find_package(OpenMP)
  endif()
else()
  message("Do not use OpenMP")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

feature_summary(WHAT ALL)

add_library(goss)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/goss/version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/goss/version.h)

set(HEADERS_goss
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/AdaptiveExplicitSolver.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/AdaptiveImplicitSolver.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/BasicImplicitEuler.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/DoubleVector.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ESDIRK23a.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ESDIRK4O32.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ExplicitEuler.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/GRL1.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/GRL2.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ImplicitEuler.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ImplicitODESolver.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ODE.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ODESolver.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ODESystemSolver.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ParameterizedODE.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/RK2.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/RK4.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/RKF32.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/RL1.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/RL2.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/constants.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/goss.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/types.h
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/version.h)

target_sources(goss PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/AdaptiveExplicitSolver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/AdaptiveImplicitSolver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/BasicImplicitEuler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ESDIRK23a.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ESDIRK4O32.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ExplicitEuler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/GRL1.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/GRL2.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ImplicitEuler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ImplicitODESolver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ODE.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ODESystemSolver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/RK2.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/RK4.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/RKF32.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/RL1.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/RL2.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/goss/ThetaSolver.cpp)

# Configure the library
set_target_properties(goss PROPERTIES PRIVATE_HEADER "${HEADERS_goss}")
target_include_directories(goss PUBLIC $<INSTALL_INTERFACE:include> "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_SOURCE_DIR}>")

if(USE_OPENMP)
  if(OpenMP_FOUND)
    target_link_libraries(goss PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()

# Set compiler flags
list(APPEND GOSS_DEVELOPER_FLAGS -O2;-g;-pipe)
list(APPEND goss_compiler_flags -Wall;-Werror;-Wextra;-Wno-comment;-pedantic)
target_compile_options(goss PRIVATE "$<$<OR:$<CONFIG:Debug>,$<CONFIG:Developer>>:${goss_compiler_flags}>")
target_compile_options(goss PRIVATE $<$<CONFIG:Developer>:${GOSS_DEVELOPER_FLAGS}>)

# Set debug definitions (private)
target_compile_definitions(goss PRIVATE $<$<OR:$<CONFIG:Debug>,$<CONFIG:Developer>>:DEBUG>)

# Install the Goss library
install(TARGETS goss
  EXPORT GossTargets
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/goss
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT RuntimeExecutables
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

# Configure CMake helpers
include(CMakePackageConfigHelpers)
write_basic_package_version_file(GossConfigVersion.cmake VERSION ${PACKAGE_VERSION}
  COMPATIBILITY AnyNewerVersion)
configure_package_config_file(GossConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/GossConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/goss)

# Install CMake files
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GossConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/GossConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/goss COMPONENT Development)
install(EXPORT GossTargets FILE GossTargets.cmake NAMESPACE Goss::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/goss)

if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(${CMAKE_SOURCE_DIR}/tests)
endif()
