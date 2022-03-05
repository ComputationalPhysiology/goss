########################################################################
# CMake build script for Google Test.
#
# To run the tests for Google Test itself on Linux, use 'make test' or
# ctest.  You can select which tests to run using 'ctest -R regex'.
# For more options, run 'ctest --help'.

include(FetchContent)
FetchContent_Declare(gtest
  QUIET
  URL https://github.com/google/googletest/archive/release-1.11.0.tar.gz
)
# configure build of googletest
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(gtest)


file(GLOB HEADERS *.h)

add_executable(testode ode.cpp)
add_executable(testodesystem odesystem.cpp)

target_link_libraries(testode
    gtest_main
    goss Boost::system Boost::thread

)

target_link_libraries(testodesystem
    gtest_main
    goss Boost::system Boost::thread

)

add_test(NAME testode
         COMMAND testode)

add_test(NAME testodesystem
         COMMAND testodesystem)
