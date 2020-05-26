if (FRACKIT_ENABLE_PYTHONBINDINGS)
    add_subdirectory(frackit)

    # modify contents (version etc) in setup file and symlink
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in
                   ${CMAKE_CURRENT_SOURCE_DIR}/setup.py)
    frackit_symlink_or_copy(FILES setup.py)

    # define target for python library
    add_custom_target(install_python
                      COMMAND "${PYTHON_EXECUTABLE}" -m pip install "${CMAKE_CURRENT_BINARY_DIR}"
                      COMMENT "Installing the python package at ${CMAKE_CURRENT_BINARY_DIR}"
                      )

    # installation rule
    install(CODE "message(\"Installing the python package wheel at ${CMAKE_CURRENT_BINARY_DIR}\")
                  execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip wheel ${CMAKE_CURRENT_BINARY_DIR}
                                  ERROR_MESSAGE \"Error installing wheel for python package at ${CMAKE_CURRENT_BINARY_DIR}\")"
           )
endif()
