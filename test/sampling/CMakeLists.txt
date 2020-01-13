###################################
# test random geometry generation #
###################################

frackit_add_test(NAME test_pointsampling SOURCES test_pointsampling.cc LABELS sampling)
frackit_add_test(NAME test_disksampling SOURCES test_disksampling.cc LABELS sampling)

set(CMAKE_BUILD_TYPE Debug)
