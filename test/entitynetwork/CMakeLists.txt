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
frackit_add_test(NAME test_generate_disk_network_cylinder
                 SOURCES test_generate_disk_network_cylinder.cc
                 LABELS entitynetwork sampling)

# test disk network creation (in shape from file)
frackit_add_test(NAME test_generate_disk_network_shape
                 SOURCES test_generate_disk_network_shape.cc
                 LABELS entitynetwork sampling
                 COMPILE_DEFINITIONS BREPFILE="${CMAKE_SOURCE_DIR}/test/entitynetwork/box.brep"
                 COMMAND ./test_generate_disk_network_shape)

set(CMAKE_BUILD_TYPE Debug)
