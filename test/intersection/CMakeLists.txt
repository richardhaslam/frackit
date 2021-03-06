###################################################
# tests for intersections of different geometries #
###################################################

# bounding box intersection
frackit_add_test(NAME test_bbox_intersection
                 SOURCES test_bbox_intersection.cc
                 LABELS intersection)

# intersection angles
frackit_add_test(NAME test_intersectionangle
                 SOURCES test_intersectionangle.cc
                 LABELS intersection)

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

# polygon - polygon
frackit_add_test(NAME test_polygon_polygon_intersection
                 SOURCES test_polygon_polygon_intersection.cc
                 LABELS intersection)

# polygon - polygon
frackit_add_test(NAME test_polygon_face_intersection
                 SOURCES test_polygon_face_intersection.cc
                 LABELS intersection)

# polygon - disk
frackit_add_test(NAME test_polygon_disk_intersection
                 SOURCES test_polygon_disk_intersection.cc
                 LABELS intersection)

# cylinder surface - disk
frackit_add_test(NAME test_cylindersurface_disk_intersection
                 SOURCES test_cylindersurface_disk_intersection.cc
                 LABELS intersection)

# cylinder surface - polygon (tests quadrilateral as well)
frackit_add_test(NAME test_cylindersurface_polygon_intersection
                 SOURCES test_cylindersurface_polygon_intersection.cc
                 LABELS intersection)

# face - disk
frackit_add_test(NAME test_disk_face_intersection
                 SOURCES test_disk_face_intersection.cc
                 LABELS intersection)

# face - quadrilateral
frackit_add_test(NAME test_quad_face_intersection
                 SOURCES test_quad_face_intersection.cc
                 LABELS intersection)

# face - face
frackit_add_test(NAME test_face_face_intersection
                 SOURCES test_face_face_intersection.cc
                 LABELS intersection)

# cpu time measurement
frackit_add_test(NAME test_intersectiontimemeasure
                 SOURCES test_intersectiontimemeasure.cc
                 LABELS intersection)

set(CMAKE_BUILD_TYPE Release)
