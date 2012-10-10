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
# Define libraries
add_library(goss ${GOSS_H} ${HEADERS} ${SOURCES})
set_target_properties(goss PROPERTIES ${GOSS_LIBRARY_PROPERTIES})

# Strip leading and trailing whitespaces
string(STRIP GOSS_TARGET_LINK_LIBRARIES "${GOSS_TARGET_LINK_LIBRARIES}")

# Add GOSS target libraries
target_link_libraries(goss ${GOSS_TARGET_LINK_LIBRARIES})
