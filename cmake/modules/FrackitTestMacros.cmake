# enable the testing suite on the CMake side.
enable_testing()

# Introduce a target that triggers the building of all tests
add_custom_target(build_tests)

# processes labels passed to frackit_add_test
function(frackit_declare_test_label)
  include(CMakeParseArguments)
  set(OPTIONS)
  set(SINGLEARGS)
  set(MULTIARGS LABELS)
  cmake_parse_arguments(arg "${OPTIONS}" "${SINGLEARGS}" "${MULTIARGS}" ${ARGN})

  if( (DEFINED arg_UNPARSED_ARGUMENTS) AND NOT ( arg_UNPARSED_ARGUMENTS STREQUAL "" ) )
    message(FATAL_ERROR "Unhandled extra arguments given to frackit_declare_test_label(): "
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
    set(target "build_${label}_tests")
    if(NOT TARGET "${target}")
      add_custom_target("${target}")
    endif()
  endforeach()
endfunction(frackit_declare_test_label)

# function to add a test
function(frackit_add_test)
  include(CMakeParseArguments)
  set(SINGLEARGS NAME TARGET TIMEOUT)
  set(MULTIARGS SOURCES COMPILE_DEFINITIONS COMPILE_FLAGS CMD_ARGS COMMAND CMAKE_GUARD LABELS)
  cmake_parse_arguments(ADDTEST "${OPTIONS}" "${SINGLEARGS}" "${MULTIARGS}" ${ARGN})

  # Check whether the parser produced any errors
  if(ADDTEST_UNPARSED_ARGUMENTS)
    message(WARNING "Unrecognized arguments ('${ADDTEST_UNPARSED_ARGUMENTS}') for frackit_add_test!")
  endif()

  # Check input for validity and apply defaults
  if(NOT ADDTEST_SOURCES AND NOT ADDTEST_TARGET)
    message(FATAL_ERROR "You need to specify either the SOURCES or the TARGET option for frackit_add_test!")
  endif()
  if(ADDTEST_SOURCES AND ADDTEST_TARGET)
    message(FATAL_ERROR "You cannot specify both SOURCES and TARGET for frackit_add_test")
  endif()
  if(NOT ADDTEST_NAME)
    message(FATAL_ERROR "You have to define a name for the test!")
  endif()
  if(NOT ADDTEST_COMMAND)
    set(ADDTEST_COMMAND ${ADDTEST_NAME})
  endif()
  if(NOT ADDTEST_TIMEOUT)
    set(ADDTEST_TIMEOUT 300)
  endif()

  # Add the executable if it is not already present
  if(ADDTEST_SOURCES)
    add_executable(${ADDTEST_NAME} ${ADDTEST_SOURCES})
    target_compile_definitions(${ADDTEST_NAME} PUBLIC ${ADDTEST_COMPILE_DEFINITIONS})
    target_compile_options(${ADDTEST_NAME} PUBLIC ${ADDTEST_COMPILE_FLAGS})

    # link GMSH library to each test
    target_link_libraries(${ADDTEST_NAME} ${ADDTEST_LINK_LIBRARIES} ${OCC_LIBS})
    set(ADDTEST_TARGET ${ADDTEST_NAME})
  endif()

  # Now add the actual test
  add_test(NAME ${ADDTEST_NAME} COMMAND "${ADDTEST_COMMAND}" ${ADDTEST_CMD_ARGS})

  # Make sure to exclude the target from all
  set_property(TARGET ${ADDTEST_TARGET} PROPERTY EXCLUDE_FROM_ALL 1)

  # make sure each label exists and its name is acceptable
  frackit_declare_test_label(LABELS ${ADDTEST_LABELS})

  # Have build_tests and build_${label}_tests depend on the given target in
  # order to trigger the build correctly
  add_dependencies(build_tests ${ADDTEST_TARGET})
  foreach(label IN LISTS ADDTEST_LABELS)
    add_dependencies(build_${label}_tests ${ADDTEST_TARGET})
  endforeach()
endfunction(frackit_add_test)
