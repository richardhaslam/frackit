#############################
# test distance computation #
#############################
frackit_add_test(NAME test_distance_linear SOURCES test_distance_linear.cc LABELS geometry distance)
frackit_add_test(NAME test_distance_disk SOURCES test_distance_disk.cc LABELS geometry distance)
frackit_add_test(NAME test_distancetoboundary SOURCES test_distancetoboundary.cc LABELS geometry distance)

set(CMAKE_BUILD_TYPE Release)
