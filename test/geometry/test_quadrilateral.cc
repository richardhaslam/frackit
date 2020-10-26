#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

#include <frackit/geometry/quadrilateral.hh>
#include <frackit/precision/precision.hh>

//! test some functionality of quadrilaterals
int main()
{
    using ctype = double;

    using Quadrilateral = Frackit::Quadrilateral<ctype, 3>;
    using Point = typename Quadrilateral::Point;

    std::vector<ctype> scales({1e-5, 1, 1e5});
    for (auto f : scales)
    {
        Quadrilateral quad(Point(0.0, 0.0, 0.0),
                           Point(f,   0.0, 0.0),
                           Point(0.0, f,   0.0),
                           Point(f,   f,   0.0));
        using std::abs;
        using std::sqrt;
        const auto eps = Frackit::Precision<ctype>::confusion()*f;
        if ( abs(quad.edge(0).length() - f) > eps )
            throw std::runtime_error("Unexpected length of edge 0");
        if ( abs(quad.edge(1).length() - f) > eps )
            throw std::runtime_error("Unexpected length of edge 1");
        if ( abs(quad.edge(2).length() - f) > eps )
            throw std::runtime_error("Unexpected length of edge 2");
        if ( abs(quad.edge(3).length() - f) > eps )
            throw std::runtime_error("Unexpected length of edge 2");
        if ( abs(quad.area() - f*f) > eps*eps )
            throw std::runtime_error("Unexpected quad area");

        // check center point
        if (!quad.center().isEqual(Point(0.5*f, 0.5*f, 0.0), eps))
            throw std::runtime_error("Center point wrong!");

        // check if quad is convex
        if (!quad.isConvex())
            throw std::runtime_error("Convexity not correct");

        // check if points are on supporting plane
        const auto plane = quad.supportingPlane();
        if (!plane.contains(quad.corner(0), eps))
            throw std::runtime_error("Point 0 not on support plane");
        if (!plane.contains(quad.corner(1), eps))
            throw std::runtime_error("Point 1 not on support plane");
        if (!plane.contains(quad.corner(2), eps))
            throw std::runtime_error("Point 2 not on support plane");
        if (!plane.contains(quad.corner(3), eps))
            throw std::runtime_error("Point 3 not on support plane");
        if (plane.contains(Point(0.0, 0.0, 1e-6*f), eps))
            throw std::runtime_error("Point should not be on plane");

        // check contains() query
        if (!quad.contains(Point(0.5*f, 0.0, 0.0), eps))
            throw std::runtime_error("contains() query 1 failed");
        if (!quad.contains(Point(0.5*f, 1e-6*f, 0.0), eps))
            throw std::runtime_error("contains() query 2 failed");
        if (quad.contains(Point(0.5*f, -1e-6*f, 0.0), eps))
            throw std::runtime_error("contains() query 3 failed");

        // make and check a non-convex quad
        Quadrilateral concave(Point(0.0, 0.0, 0.0),
                              Point(f,   0.0, 0.0),
                              Point(0.0, f,   0.0),
                              Point(0.5*f, 0.25*f,   0.0));

        if (concave.isConvex())
            throw std::runtime_error("Concave quad not detected");

        // check a quad that is degenerated to a triangle
        Quadrilateral tria(Point(0.0, 0.0, 0.0),
                           Point(f,   0.0, 0.0),
                           Point(0.0, f,   0.0),
                           Point(0.5*f, 0.5*f,   0.0));

        if (!tria.isConvex())
            throw std::runtime_error("Triangle detected as non-convex");

        // check a quad that is slightly non-convex
        Quadrilateral slightlyConcave(Point(0.0, 0.0, 0.0),
                                      Point(f,   0.0, 0.0),
                                      Point(0.0, f,   0.0),
                                      Point(0.5*f, 0.5*f - 1e-6*f,   0.0));

        if (slightlyConcave.isConvex())
            throw std::runtime_error("Slightly concave quad not detected");

        // check a quad that is slightly non-convex
        Quadrilateral slightlyConvex(Point(0.0, 0.0, 0.0),
                                     Point(f,   0.0, 0.0),
                                     Point(0.0, f,   0.0),
                                     Point(0.5*f, 0.5*f + 1e-6*f,   0.0));

        if (!slightlyConvex.isConvex())
            throw std::runtime_error("Slightly convex quad not detected");
    }

    std::cout << "All tests passed" << std::endl;
    return 0;
}
