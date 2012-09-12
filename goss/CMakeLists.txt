set(GOSS_H goss.h)

#------------------------------------------------------------------------------
# GOSS source directories

# Headers in this directory
file(GLOB HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h)

# Sources in this directory
file(GLOB SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

#------------------------------------------------------------------------------
# Install header files

install(FILES ${GOSS_H} DESTINATION ${GOSS_INCLUDE_DIR} COMPONENT Development)

foreach(_HEADER ${HEADERS})
  list(FIND ${GOSS_H} ${_HEADER} INDEX_EXCLUDES)
  if (${INDEX_EXCLUDES} LESS 0)
    install(FILES ${_HEADER} DESTINATION ${GOSS_INCLUDE_DIR}/goss COMPONENT Development)
  endif()
endforeach()

#------------------------------------------------------------------------------
# Set compiler flags, include directories and library dependencies

# Add compiler include directories
# include_directories(${GOSS_SOURCE_DIR} ${GOSS_DEP_INCLUDE_DIRECTORIES})
# include_directories(SYSTEM ${GOSS_DEP_SYSTEM_INCLUDE_DIRECTORIES})

# Add CXX defintions
add_definitions(${GOSS_CXX_DEFINITIONS})
add_definitions(-DGOSS_VERSION="${GOSS_VERSION}")

# Add flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GOSS_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GOSS_LINK_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${GOSS_LINK_FLAGS}")

# Define libraries
add_library(goss ${GOSS_H} ${HEADERS} ${SOURCES})
set_target_properties(goss PROPERTIES ${GOSS_LIBRARY_PROPERTIES})

# Strip leading and trailing whitespaces
string(STRIP GOSS_TARGET_LINK_LIBRARIES "${GOSS_TARGET_LINK_LIBRARIES}")

# Add GOSS target libraries
target_link_libraries(goss ${GOSS_TARGET_LINK_LIBRARIES})

# Convert libraries to -L<libdir> -l<lib> form
foreach(_lib ${GOSS_TARGET_LINK_LIBRARIES})
  # Add -Wl,option directives
  if ("${_lib}" MATCHES "-Wl,[^ ]*")
     set(PKG_LINKFLAGS "${_lib} ${PKG_LINKFLAGS}")
  else()
    string(REGEX REPLACE "(.?:?/[^ ]*)/lib([^ ]*)\\.(a|so|dylib|dll)" "-L\\1 -l\\2"
      _linkflags
      "${_lib}"
      )

    # Add libraries that matches the form -L<libdir> -l<lib>
    if ("${_linkflags}" MATCHES "-L.+ -l.+")
        set(PKG_LINKFLAGS "${_linkflags} ${PKG_LINKFLAGS}")
    endif()
  endif()
endforeach()

message(STATUS "LINKFLAGS ${PKG_LINKFLAGS}")

# Remove duplicated link flags
if (PKG_LINKFLAGS)
  separate_arguments(PKG_LINKFLAGS)
  list(REMOVE_DUPLICATES PKG_LINKFLAGS)
  string(REPLACE ";" " " PKG_LINKFLAGS "${PKG_LINKFLAGS}")
endif()

# Add additional link flags
foreach(_linkflag ${GOSS_LINK_FLAGS})
  set(PKG_LINKFLAGS "${PKG_LINKFLAGS} ${_linkflag}")
endforeach()
