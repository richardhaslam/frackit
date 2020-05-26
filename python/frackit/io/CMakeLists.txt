frackit_symlink_or_copy(FILES __init__.py)
pybind11_add_module(_io _io.cc)
target_link_libraries(_io PRIVATE ${OCC_LIBS})
