cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(goss VERSION 0.1.0)

# Include dir
include_directories(/usr/local/include)

# Src
AUX_SOURCE_DIRECTORY(src SRC_FILES)

# Headers
set(PROJECT_SOURCE_DIR "src")
set(PROJECT_INCLUDE_DIR "include/goss_bits")

# Source files
set(SOURCE_FILES
    ${PROJECT_INCLUDE_DIR}/AdaptiveExplicitSolver.h
    ${PROJECT_INCLUDE_DIR}/AdaptiveImplicitSolver.h
    ${PROJECT_INCLUDE_DIR}/BasicImplicitEuler.h
    ${PROJECT_INCLUDE_DIR}/DoubleVector.h
    ${PROJECT_INCLUDE_DIR}/ESDIRK23a.h
    ${PROJECT_INCLUDE_DIR}/ESDIRK4O32.h
    ${PROJECT_INCLUDE_DIR}/ExplicitEuler.h
    ${PROJECT_INCLUDE_DIR}/GRL1.h
    ${PROJECT_INCLUDE_DIR}/GRL2.h
    ${PROJECT_INCLUDE_DIR}/ImplicitEuler.h
    ${PROJECT_INCLUDE_DIR}/ImplicitODESolver.h
    ${PROJECT_INCLUDE_DIR}/LogLevel.h
    ${PROJECT_INCLUDE_DIR}/LogManager.h
    ${PROJECT_INCLUDE_DIR}/LogStream.h
    ${PROJECT_INCLUDE_DIR}/Logger.h
    ${PROJECT_INCLUDE_DIR}/ODE.h
    ${PROJECT_INCLUDE_DIR}/ODESolver.h
    ${PROJECT_INCLUDE_DIR}/ODESystemSolver.h
    ${PROJECT_INCLUDE_DIR}/Parameter.h
    ${PROJECT_INCLUDE_DIR}/ParameterizedODE.h
    ${PROJECT_INCLUDE_DIR}/Parameters.h
    ${PROJECT_INCLUDE_DIR}/Progress.h
    ${PROJECT_INCLUDE_DIR}/RK2.h
    ${PROJECT_INCLUDE_DIR}/RK4.h
    ${PROJECT_INCLUDE_DIR}/RKF32.h
    ${PROJECT_INCLUDE_DIR}/RL1.h
    ${PROJECT_INCLUDE_DIR}/RL2.h
    ${PROJECT_INCLUDE_DIR}/Table.h
    ${PROJECT_INCLUDE_DIR}/ThetaSolver.h
    ${PROJECT_INCLUDE_DIR}/Timer.h
    ${PROJECT_INCLUDE_DIR}/constants.h
    ${PROJECT_INCLUDE_DIR}/goss.h
    ${PROJECT_INCLUDE_DIR}/log.h
    ${PROJECT_INCLUDE_DIR}/timing.h
    ${PROJECT_INCLUDE_DIR}/types.h
    ${PROJECT_SOURCE_DIR}/AdaptiveExplicitSolver.cpp
    ${PROJECT_SOURCE_DIR}/AdaptiveImplicitSolver.cpp
    ${PROJECT_SOURCE_DIR}/BasicImplicitEuler.cpp
    ${PROJECT_SOURCE_DIR}/ESDIRK23a.cpp
    ${PROJECT_SOURCE_DIR}/ESDIRK4O32.cpp
    ${PROJECT_SOURCE_DIR}/ExplicitEuler.cpp
    ${PROJECT_SOURCE_DIR}/GRL1.cpp
    ${PROJECT_SOURCE_DIR}/GRL2.cpp
    ${PROJECT_SOURCE_DIR}/ImplicitEuler.cpp
    ${PROJECT_SOURCE_DIR}/ImplicitODESolver.cpp
    ${PROJECT_SOURCE_DIR}/LogManager.cpp
    ${PROJECT_SOURCE_DIR}/LogStream.cpp
    ${PROJECT_SOURCE_DIR}/Logger.cpp
    ${PROJECT_SOURCE_DIR}/ODE.cpp
    ${PROJECT_SOURCE_DIR}/ODESystemSolver.cpp
    ${PROJECT_SOURCE_DIR}/Parameter.cpp
    ${PROJECT_SOURCE_DIR}/Parameters.cpp
    ${PROJECT_SOURCE_DIR}/Progress.cpp
    ${PROJECT_SOURCE_DIR}/RK2.cpp
    ${PROJECT_SOURCE_DIR}/RK4.cpp
    ${PROJECT_SOURCE_DIR}/RKF32.cpp
    ${PROJECT_SOURCE_DIR}/RL1.cpp
    ${PROJECT_SOURCE_DIR}/RL2.cpp
    ${PROJECT_SOURCE_DIR}/Table.cpp
    ${PROJECT_SOURCE_DIR}/ThetaSolver.cpp
    ${PROJECT_SOURCE_DIR}/log.cpp
    ${PROJECT_SOURCE_DIR}/timing.cpp
)

# Set up such that XCode organizes the files correctly
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES})

# Add library
add_library(goss SHARED ${SOURCE_FILES})

# Include directories
target_include_directories(goss PRIVATE include/goss_bits)

TARGET_LINK_LIBRARIES(goss LINK_PUBLIC ${Boost_LIBRARIES})

# Install
install(TARGETS goss DESTINATION lib)

# Install the headers
install(FILES include/goss_bits DESTINATION include)

# Create base directory
install(DIRECTORY include/goss_bits DESTINATION include)