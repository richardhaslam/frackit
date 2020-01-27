frackit_add_application(NAME example3
                        SOURCES example3.cc
                        COMPILE_DEFINITIONS BREPFILE="${CMAKE_SOURCE_DIR}/appl/example3/layers.brep"
                        LABELS examples)

set(CMAKE_BUILD_TYPE Debug)
