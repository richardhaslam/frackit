###################################
# test random geometry generation #
###################################

frackit_add_test(NAME test_status SOURCES test_status.cc LABELS sampling)
frackit_add_test(NAME test_pointsampling SOURCES test_pointsampling.cc LABELS sampling)
frackit_add_test(NAME test_disksampling SOURCES test_disksampling.cc LABELS sampling)
frackit_add_test(NAME test_quadrilateralsampling SOURCES test_quadrilateralsampling.cc LABELS sampling)
frackit_add_test(NAME test_polygonsampling SOURCES test_polygonsampling.cc LABELS sampling)

set(CMAKE_BUILD_TYPE Release)
