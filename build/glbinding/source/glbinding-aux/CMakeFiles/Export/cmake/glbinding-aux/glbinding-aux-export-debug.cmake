#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "glbinding::glbinding-aux" for configuration "Debug"
set_property(TARGET glbinding::glbinding-aux APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(glbinding::glbinding-aux PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/glbinding-auxd.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS glbinding::glbinding-aux )
list(APPEND _IMPORT_CHECK_FILES_FOR_glbinding::glbinding-aux "${_IMPORT_PREFIX}/lib/glbinding-auxd.lib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
