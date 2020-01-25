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

# quadrilateral - line
frackit_add_test(NAME test_quad_line_intersection
                 SOURCES test_quad_line_intersection.cc
                 LABELS intersection)

# disk - disk
frackit_add_test(NAME test_disk_disk_intersection
                 SOURCES test_disk_disk_intersection.cc
                 LABELS intersection)

# quadrilateral - quadrilateral
frackit_add_test(NAME test_quad_quad_intersection
                 SOURCES test_quad_quad_intersection.cc
                 LABELS intersection)

# quadrilateral - disk
frackit_add_test(NAME test_quad_disk_intersection
                 SOURCES test_quad_disk_intersection.cc
                 LABELS intersection)

# cylinder surface - disk
frackit_add_test(NAME test_cylindersurface_disk_intersection
                 SOURCES test_cylindersurface_disk_intersection.cc
                 LABELS intersection)

# cylinder surface - quadrilateral
frackit_add_test(NAME test_cylindersurface_quad_intersection
                 SOURCES test_cylindersurface_quad_intersection.cc
                 LABELS intersection)

# face - disk
frackit_add_test(NAME test_disk_face_intersection
                 SOURCES test_disk_face_intersection.cc
                 LABELS intersection)

# face - quadrilateral
frackit_add_test(NAME test_quad_face_intersection
                 SOURCES test_quad_face_intersection.cc
                 LABELS intersection)

# shell - disk
frackit_add_test(NAME test_disk_shell_intersection
                 SOURCES test_disk_shell_intersection.cc
                 LABELS intersection)

set(CMAKE_BUILD_TYPE Debug)
