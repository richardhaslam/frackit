######################################
# test functionalities of geometries #
######################################
frackit_add_test(NAME test_point SOURCES test_point.cc LABELS geometry)
frackit_add_test(NAME test_direction SOURCES test_direction.cc LABELS geometry)
frackit_add_test(NAME test_plane SOURCES test_plane.cc LABELS geometry)
frackit_add_test(NAME test_circle SOURCES test_circle.cc LABELS geometry)
frackit_add_test(NAME test_cylindricalsurface SOURCES test_cylindricalsurface.cc LABELS geometry)
frackit_add_test(NAME test_cylinder SOURCES test_cylinder.cc LABELS geometry)
frackit_add_test(NAME test_hollowcylinder SOURCES test_hollowcylinder.cc LABELS geometry)
frackit_add_test(NAME test_disk SOURCES test_disk.cc LABELS geometry)
frackit_add_test(NAME test_ellipsearc SOURCES test_ellipsearc.cc LABELS geometry)
frackit_add_test(NAME test_triangle SOURCES test_triangle.cc LABELS geometry)
frackit_add_test(NAME test_quadrilateral SOURCES test_quadrilateral.cc LABELS geometry)
frackit_add_test(NAME test_box SOURCES test_box.cc LABELS geometry)
frackit_add_test(NAME test_name SOURCES test_name.cc LABELS geometry)

set(CMAKE_BUILD_TYPE Release)
