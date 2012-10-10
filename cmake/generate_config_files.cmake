#------------------------------------------------------------------------------
# Generate pkg-config file and install it

# Define packages that should be required by pkg-config file
set(PKG_REQUIRES "")

# Convert include dirs to -I<incdir> form
foreach(_inc_dir ${GOSS_DEP_INCLUDE_DIRECTORIES})
  set(PKG_INCLUDES "-I${_inc_dir} ${PKG_INCLUDES}")
endforeach()

foreach(_inc_dir ${GOSS_DEP_SYSTEM_INCLUDE_DIRECTORIES})
  set(PKG_INCLUDES "-I${_inc_dir} ${PKG_INCLUDES}")
endforeach()

# Convert compiler flags and definitions into space separated strings
string(REPLACE ";" " " PKG_CXXFLAGS "${CMAKE_CXX_FLAGS}")
string(REPLACE ";" " " PKG_LINKFLAGS "${CMAKE_EXE_LINKER_FLAGS}")
string(REPLACE ";" " " PKG_DEFINITIONS "${GOSS_CXX_DEFINITIONS}")

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

# Remove duplicated link flags
separate_arguments(PKG_LINKFLAGS)
list(REMOVE_DUPLICATES PKG_LINKFLAGS)
string(REPLACE ";" " " PKG_LINKFLAGS "${PKG_LINKFLAGS}")

# Add additional link flags
foreach(_linkflag ${GOSS_LINK_FLAGS})
  set(PKG_LINKFLAGS "${PKG_LINKFLAGS} ${_linkflag}")
endforeach()

# Configure and install pkg-config file
configure_file(${GOSS_CMAKE_DIR}/templates/goss.pc.in ${CMAKE_BINARY_DIR}/goss.pc @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/goss.pc
  DESTINATION ${GOSS_PKGCONFIG_DIR}
  COMPONENT Development
  )

#------------------------------------------------------------------------------
# Generate CMake config files (goss-config{,-version}.cmake)

# Set library location
get_target_property(GOSS_LIBRARY_LOCATION goss LOCATION)
get_filename_component(GOSS_LIBRARY_FILENAME ${GOSS_LIBRARY_LOCATION} NAME)
set(GOSS_LIBRARY "${CMAKE_INSTALL_PREFIX}/${GOSS_LIB_DIR}/${GOSS_LIBRARY_FILENAME}")

configure_file(${GOSS_CMAKE_DIR}/templates/goss-config.cmake.in
  ${CMAKE_BINARY_DIR}/goss/goss-config.cmake @ONLY)
configure_file(${GOSS_CMAKE_DIR}/templates/goss-config-version.cmake.in
  ${CMAKE_BINARY_DIR}/goss/goss-config-version.cmake @ONLY)
install(
  FILES
    ${CMAKE_BINARY_DIR}/goss/goss-config.cmake
    ${CMAKE_BINARY_DIR}/goss/goss-config-version.cmake
  DESTINATION ${GOSS_SHARE_DIR}/cmake
  COMPONENT Development
  )

#------------------------------------------------------------------------------
# Generate and install helper file goss.conf

# FIXME: Can CMake provide the library path name variable?
if (APPLE)
  set(OS_LIBRARY_PATH_NAME "DYLD_LIBRARY_PATH")
else()
  set(OS_LIBRARY_PATH_NAME "LD_LIBRARY_PATH")
endif()

# FIXME: not cross-platform compatible
# Create and install goss.conf file
configure_file(${GOSS_CMAKE_DIR}/templates/goss.conf.in
               ${CMAKE_BINARY_DIR}/goss.conf @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/goss.conf
        DESTINATION ${GOSS_SHARE_DIR}
        COMPONENT Development)

#------------------------------------------------------------------------------
# Generate and install helper file goss-version

# FIXME: not cross-platform compatible
# Create and install goss-version file
configure_file(${GOSS_CMAKE_DIR}/templates/goss-version.in
               ${CMAKE_BINARY_DIR}/goss-version @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/goss-version
        DESTINATION ${GOSS_BIN_DIR}
	PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE
        COMPONENT RuntimeExecutables)

