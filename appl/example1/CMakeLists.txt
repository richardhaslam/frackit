frackit_add_application(NAME example1
                        SOURCES example1.cc
                        LABELS example)

frackit_symlink_or_copy(FILES example1.py)

set(CMAKE_BUILD_TYPE Release)
