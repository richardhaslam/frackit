###################################################
# tests for intersections of different geometries #
###################################################

# segment - segment
frackit_add_test(NAME test_segment_segment_intersection
                 SOURCES test_segment_segment_intersection.cc
                 LABELS intersection)

# plane - plane
frackit_add_test(NAME test_plane_plane_intersection
                 SOURCES test_plane_plane_intersection.cc
                 LABELS intersection)

# disk - line
frackit_add_test(NAME test_disk_line_intersection
                 SOURCES test_disk_line_intersection.cc
                 LABELS intersection)

# disk - disk
frackit_add_test(NAME test_disk_disk_intersection
                 SOURCES test_disk_disk_intersection.cc
                 LABELS intersection)

# cylinder surface - disk
frackit_add_test(NAME test_cylindersurface_disk_intersection
                 SOURCES test_cylindersurface_disk_intersection.cc
                 LABELS intersection)

set(CMAKE_BUILD_TYPE Debug)
