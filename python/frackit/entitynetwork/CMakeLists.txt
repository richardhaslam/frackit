frackit_symlink_or_copy(FILES __init__.py)
pybind11_add_module(_entitynetwork _entitynetwork.cc)
target_link_libraries(_entitynetwork PRIVATE ${OCC_LIBS})
