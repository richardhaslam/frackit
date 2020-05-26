frackit_symlink_or_copy(FILES __init__.py)
pybind11_add_module(_occutilities _occutilities.cc)
target_link_libraries(_occutilities PRIVATE ${OCC_LIBS})
