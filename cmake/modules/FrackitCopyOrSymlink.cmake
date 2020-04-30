
# function to create a symbolic link between files in the source and build directories
function(frackit_symlink_or_copy)

  # skip this if source and binary dir are the same
  if( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR} )
    return()
  endif()

  # parse arguments
  include(CMakeParseArguments)
  cmake_parse_arguments(ARG "" "DESTINATION" "FILES" ${ARGN})
  if (ARG_UNPARSED_ARGUMENTS)
    message(WARNING "There were issues parsing your arguments to frackit_symlink_or_copy_to_source)!")
  endif ()

  # handle all given files
  foreach (f ${ARG_FILES})
    # check whether a speficic destination is desired
    if (ARG_DESTINATION)
      set(destination "${ARG_DESTINATION}/")
    else ()
      set(destination "")
    endif ()

    # On windows, copy instead of linking
    if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
      if(NOT DEFINED FRACKIT_SYMLINK_WARNING)
        message(WARNING "Symlinks are not supported on your platform! Doing copies instead.")
        set(FRACKIT_SYMLINK_WARNING)
      endif()
      execute_process(COMMAND ${CMAKE_COMMAND} "-E" "copy" "${CMAKE_CURRENT_SOURCE_DIR}/${f}" "${CMAKE_CURRENT_BINARY_DIR}/${destination}${f}")
    else ()
      # create symlink
      execute_process(COMMAND ${CMAKE_COMMAND} "-E" "create_symlink" "${CMAKE_CURRENT_SOURCE_DIR}/${f}" "${CMAKE_CURRENT_BINARY_DIR}/${destination}${f}")
    endif ()
  endforeach ()

endfunction(frackit_symlink_or_copy)
