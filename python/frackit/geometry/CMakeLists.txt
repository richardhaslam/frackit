frackit_symlink_or_copy(FILES __init__.py)
pybind11_add_module(_geometry _geometry.cc)
target_link_libraries(_geometry PRIVATE ${OCC_LIBS})
