frackit_symlink_or_copy(FILES __init__.py)
pybind11_add_module(_intersection _intersection.cc)
target_link_libraries(_intersection PRIVATE ${OCC_LIBS})
