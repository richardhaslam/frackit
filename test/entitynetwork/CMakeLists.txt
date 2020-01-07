##########################################
# tests functionality on entity networks #
##########################################

# test network constaints
frackit_add_test(NAME test_constraints_disk
                 SOURCES test_constraints_disk.cc
                 LABELS entitynetwork constraints)

set(CMAKE_BUILD_TYPE Debug)
