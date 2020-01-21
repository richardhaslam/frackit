# test write of networks to gmsh output format
frackit_add_test(NAME test_write_gmsh
                 SOURCES test_write_gmsh.cc
                 LABELS io gmsh)

set(CMAKE_BUILD_TYPE Debug)
