#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "glbinding::glbinding-aux" for configuration "RelWithDebInfo"
set_property(TARGET glbinding::glbinding-aux APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(glbinding::glbinding-aux PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "CXX"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/glbinding-aux.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS glbinding::glbinding-aux )
list(APPEND _IMPORT_CHECK_FILES_FOR_glbinding::glbinding-aux "${_IMPORT_PREFIX}/lib/glbinding-aux.lib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
