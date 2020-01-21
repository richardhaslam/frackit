frackit_add_test(NAME example_layers
                 SOURCES example_layers.cc
                 COMPILE_DEFINITIONS BREPFILE="${CMAKE_SOURCE_DIR}/examples/layers/layers.brep"
                 LABELS examples)

set(CMAKE_BUILD_TYPE Debug)
