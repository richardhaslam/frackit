frackit_add_application(NAME example4
                        SOURCES example4.cc
                        LABELS example)

frackit_symlink_or_copy(FILES example4.py)

set(CMAKE_BUILD_TYPE Release)
