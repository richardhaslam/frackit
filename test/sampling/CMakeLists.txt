###################################
# test random geometry generation #
###################################

# test point sampling
frackit_add_test(NAME test_pointsampling
                 SOURCES test_pointsampling.cc
                 LABELS sampling)

set(CMAKE_BUILD_TYPE Debug)
