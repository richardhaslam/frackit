frackit_symlink_or_copy(FILES __init__.py)
pybind11_add_module(_precision _precision.cc)
target_link_libraries(_precision PRIVATE ${OCC_LIBS})
