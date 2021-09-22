########################################################################
# CMake build script for Google Test.
#
# To run the tests for Google Test itself on Linux, use 'make test' or
# ctest.  You can select which tests to run using 'ctest -R regex'.
# For more options, run 'ctest --help'.

include(FetchContent)
FetchContent_Declare(gtest
  QUIET
  URL https://github.com/google/googletest/archive/release-1.10.0.tar.gz
)
# configure build of googletest
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(gtest)


file(GLOB HEADERS *.h)

add_executable(testode ode.cpp ${HEADERS})
add_executable(testodesystem odesystem.cpp ${HEADERS})

target_link_libraries(testode
    goss
    gtest_main
    Boost::boost

)

target_link_libraries(testodesystem
    goss
    gtest_main
    Boost::boost

)

add_test(NAME testode
         COMMAND testode)

add_test(NAME testodesystem
         COMMAND testodesystem)

