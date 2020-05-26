frackit_symlink_or_copy(FILES __init__.py)
pybind11_add_module(_magnitude _magnitude.cc)
target_link_libraries(_magnitude PRIVATE ${OCC_LIBS})
