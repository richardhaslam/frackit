# Introduce a target that triggers the building of all applications
add_custom_target(build_applications)

# processes labels passed to frackit_add_application
function(frackit_declare_application_label)
  include(CMakeParseArguments)
  set(OPTIONS)
  set(SINGLEARGS)
  set(MULTIARGS LABELS)
  cmake_parse_arguments(arg "${OPTIONS}" "${SINGLEARGS}" "${MULTIARGS}" ${ARGN})

  if( (DEFINED arg_UNPARSED_ARGUMENTS) AND NOT ( arg_UNPARSED_ARGUMENTS STREQUAL "" ) )
    message(FATAL_ERROR "Unhandled extra arguments given to frackit_declare_application_label(): "
      "<${arg_UNPARSED_ARGUMENTS}>")
  endif()

  foreach(label IN LISTS arg_LABELS)
    # Make sure the label is not empty, and does not contain any funny
    # characters, in particular regex characters
    if(NOT (label MATCHES "[-_0-9a-zA-Z]+"))
      message(FATAL_ERROR "Refusing to add label \"${label}\" since it is "
        "empty or contains funny characters (characters other than "
        "alphanumeric ones and \"-\" or \"_\"; the intent of this restriction "
        "is to make construction of the argument to \"ctest -L\" easier")
    endif()
    set(target "build_${label}_applications")
    if(NOT TARGET "${target}")
      add_custom_target("${target}")
    endif()
  endforeach()
endfunction(frackit_declare_application_label)

# function to add an application
function(frackit_add_application)
  include(CMakeParseArguments)
  set(SINGLEARGS NAME TARGET TIMEOUT)
  set(MULTIARGS SOURCES COMPILE_DEFINITIONS COMPILE_FLAGS CMD_ARGS COMMAND CMAKE_GUARD LABELS)
  cmake_parse_arguments(ADDAPPL "${OPTIONS}" "${SINGLEARGS}" "${MULTIARGS}" ${ARGN})

  # Check whether the parser produced any errors
  if(ADDAPPL_UNPARSED_ARGUMENTS)
    message(WARNING "Unrecognized arguments ('${ADDAPPL_UNPARSED_ARGUMENTS}') for frackit_add_application!")
  endif()

  # Check input for validity and apply defaults
  if(NOT ADDAPPL_SOURCES AND NOT ADDAPPL_TARGET)
    message(FATAL_ERROR "You need to specify either the SOURCES or the TARGET option for frackit_add_application!")
  endif()
  if(ADDAPPL_SOURCES AND ADDAPPL_TARGET)
    message(FATAL_ERROR "You cannot specify both SOURCES and TARGET for frackit_add_application")
  endif()
  if(NOT ADDAPPL_NAME)
    message(FATAL_ERROR "You have to define a name for the application!")
  endif()
  if(NOT ADDAPPL_COMMAND)
    set(ADDAPPL_COMMAND ${ADDAPPL_NAME})
  endif()
  if(NOT ADDAPPL_TIMEOUT)
    set(ADDAPPL_TIMEOUT 3600)
  endif()

  # Add the executable if it is not already present
  if(ADDAPPL_SOURCES)
    add_executable(${ADDAPPL_NAME} ${ADDAPPL_SOURCES})
    target_compile_definitions(${ADDAPPL_NAME} PUBLIC ${ADDAPPL_COMPILE_DEFINITIONS})
    target_compile_options(${ADDAPPL_NAME} PUBLIC ${ADDAPPL_COMPILE_FLAGS})

    # link OCC library to each application
    target_link_libraries(${ADDAPPL_NAME} ${ADDAPPL_LINK_LIBRARIES} ${OCC_LIBS})
    set(ADDAPPL_TARGET ${ADDAPPL_NAME})
  endif()

  # make sure each label exists and its name is acceptable
  frackit_declare_application_label(LABELS ${ADDAPPL_LABELS})

  # Make sure to exclude the target from all
  set_property(TARGET ${ADDAPPL_TARGET} PROPERTY EXCLUDE_FROM_ALL 1)


  # Have build_applications and build_${label}_applications depend on the given target in
  # order to trigger the build correctly
  add_dependencies(build_applications ${ADDAPPL_TARGET})
  foreach(label IN LISTS ADDAPPL_LABELS)
    add_dependencies(build_${label}_applications ${ADDAPPL_TARGET})
  endforeach()

endfunction(frackit_add_application)
