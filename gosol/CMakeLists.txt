set(GOSOL_H gosol.h)

#------------------------------------------------------------------------------
# GOSOL source directories

# Headers in this directory
file(GLOB HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h)

# Sources in this directory
file(GLOB SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

#------------------------------------------------------------------------------
# Install header files

install(FILES ${GOSOL_H} DESTINATION ${GOSOL_INCLUDE_DIR} COMPONENT Development)

foreach(_HEADER ${HEADERS})
  list(FIND ${GOSOL_H} ${_HEADER} INDEX_EXCLUDES)
  if (${INDEX_EXCLUDES} LESS 0)
    install(FILES ${_HEADER} DESTINATION ${GOSOL_INCLUDE_DIR}/gosol COMPONENT Development)
  endif()
endforeach()

#------------------------------------------------------------------------------
# Set compiler flags, include directories and library dependencies

# Add compiler include directories
# include_directories(${GOSOL_SOURCE_DIR} ${GOSOL_DEP_INCLUDE_DIRECTORIES})
# include_directories(SYSTEM ${GOSOL_DEP_SYSTEM_INCLUDE_DIRECTORIES})

# Add CXX defintions
add_definitions(${GOSOL_CXX_DEFINITIONS})
add_definitions(-DGOSOL_VERSION="${GOSOL_VERSION}")

# Add flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GOSOL_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GOSOL_LINK_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${GOSOL_LINK_FLAGS}")

# Define libraries
add_library(gosol ${GOSOL_H} ${HEADERS} ${SOURCES})
set_target_properties(gosol PROPERTIES ${GOSOL_LIBRARY_PROPERTIES})

# Strip leading and trailing whitespaces
string(STRIP GOSOL_TARGET_LINK_LIBRARIES "${GOSOL_TARGET_LINK_LIBRARIES}")

# Add GOSOL target libraries
target_link_libraries(gosol ${GOSOL_TARGET_LINK_LIBRARIES})

# Convert libraries to -L<libdir> -l<lib> form
foreach(_lib ${GOSOL_TARGET_LINK_LIBRARIES})
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
foreach(_linkflag ${GOSOL_LINK_FLAGS})
  set(PKG_LINKFLAGS "${PKG_LINKFLAGS} ${_linkflag}")
endforeach()
