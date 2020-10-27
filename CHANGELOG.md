Differences between Frackit 1.2 and Frackit 1.1
===============================================

### Improvements and Enhancements
- An implementation of `Sphere` has been added to the classes of available geometries
- Three-dimensional `Polygon` class has been implemented
- A sampler class for 3d polygons, `PolygonSampler`, has been added
- The `QuadrilateralSampler` class now internally uses `PolygonSampler`, leading to more variance in the output geometries

### Deprecated classes to be removed after 1.3
- Since `QuadrilateralSampler` is now a special case of the more generic `PolygonSampler`, the constructor and underlying
parameter distributions have been changed. The old constructor/behaviour is still supported, although throwing a deprecation
warning. Note that in case you were using `QuadrilateralSampler` with your own custom traits class, your code might not be compatible.
Backwards compatibility could only be achieved for usage with the default traits.

Differences between Frackit 1.1 and Frackit 1.0
===============================================

### Improvements and Enhancements

- Python bindings are now available for almost all functionality
- Python code to run all examples were added
