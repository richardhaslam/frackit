frackit_add_application(NAME example2
                        SOURCES example2.cc
                        LABELS example)

frackit_symlink_or_copy(FILES example2.py)

set(CMAKE_BUILD_TYPE Release)
