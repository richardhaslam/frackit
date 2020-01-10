##########################################
# tests functionality on entity networks #
##########################################

# test network constaints (disks)
frackit_add_test(NAME test_constraints_disk
                 SOURCES test_constraints_disk.cc
                 LABELS entitynetwork constraints)

# test network constaints (cylinder surface)
frackit_add_test(NAME test_constraints_cylindersurface
                 SOURCES test_constraints_cylindersurface.cc
                 LABELS entitynetwork constraints)

# test disk network creation (in cylinder)
frackit_add_test(NAME test_generate_disk_network
                 SOURCES test_generate_disk_network.cc
                 LABELS entitynetwork sampling)

set(CMAKE_BUILD_TYPE Debug)
