frackit_add_application(NAME example3
                        SOURCES example3.cc
                        COMPILE_DEFINITIONS BREPFILE="layers.brep"
                        LABELS example)

frackit_symlink_or_copy(FILES example3.py layers.brep)

set(CMAKE_BUILD_TYPE Release)
