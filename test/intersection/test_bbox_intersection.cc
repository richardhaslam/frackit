#include <stdexcept>

#include <frackit/geometry/boundingbox.hh>
#include <frackit/intersection/bboxintersection.hh>

//! test bounding box intersections
//! TODO: Implement and test also the boxes that result from intersection
int main()
{
    using namespace Frackit;
    using ctype = double;
    using BBox = BoundingBox<ctype>;

    std::vector<ctype> scales({1.0e-8, 1.0e-5, 1, 1e5, 1e8});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;

        const BBox unitBox(0.0, 0.0, 0.0, f, f, f);

        // identical box
        if (!doIntersect(unitBox, unitBox))
            throw std::runtime_error("Self-intersection not detected");
        std::cout << "Test 1 passed" << std::endl;

        // immersed box
        const BBox immersed(0.2*f, 0.2*f, 0.5*f, 0.8*f, 0.8*f, 0.8*f);
        if (!doIntersect(unitBox, immersed))
            throw std::runtime_error("Immersed box not detected");
        std::cout << "Test 2.1 passed" << std::endl;
        if (!doIntersect(immersed, unitBox))
            throw std::runtime_error("Immersed box not detected");
        std::cout << "Test 2.1 passed" << std::endl;

        // face touching box
        const BBox faceTouching(f, 0.2*f, 0.5*f, 2.0*f, 2.0*f, 2.0*f);
        if (!doIntersect(unitBox, faceTouching))
            throw std::runtime_error("Face touching box not detected");
        std::cout << "Test 3.1 passed" << std::endl;
        if (!doIntersect(faceTouching, unitBox))
            throw std::runtime_error("Face touching box not detected");
        std::cout << "Test 3.1 passed" << std::endl;

        // corner touching box
        const BBox cornerTouching(f, f, f, 2.0*f, 2.0*f, 2.0*f);
        if (!doIntersect(unitBox, cornerTouching))
            throw std::runtime_error("Face touching box not detected");
        std::cout << "Test 3.1 passed" << std::endl;
        if (!doIntersect(cornerTouching, unitBox))
            throw std::runtime_error("Face touching box not detected");
        std::cout << "Test 3.1 passed" << std::endl;

        std::cout << "All tests passed\n" << std::endl;
    }

    return 0;
}
