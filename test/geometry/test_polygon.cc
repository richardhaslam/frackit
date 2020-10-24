#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

#include <frackit/geometry/polygon.hh>
#include <frackit/precision/precision.hh>

//! test some functionality of polygons
int main()
{
    using ctype = double;

    using Polygon = Frackit::Polygon<ctype, 3>;
    using Point = typename Polygon::Point;

    std::vector<ctype> scales({1e-5, 1, 1e5});

    for (auto f : scales)
    {
        Polygon quad(std::vector<Point>{{Point(0.0, 0.0, 0.0),
                                         Point(f,   0.0, 0.0),
                                         Point(f, f,   0.0),
                                         Point(0.0,   f,   0.0)}});
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

        Polygon polygon(std::vector<Point>{{Point(0.0, 0.0, 0.0),
                                            Point(f,   0.0, 0.0),
                                            Point(1.5*f, 0.5*f, 0.0),
                                            Point(f,   f,   0.0),
                                            Point(0, f,   0.0),
                                            Point(-0.5*f, 0.5*f, 0.0)}});

        // check if points are on supporting plane
        const auto polyPlane = polygon.supportingPlane();
        for (unsigned int i = 0; i < polygon.numCorners(); ++i)
            if (!polyPlane.contains(polygon.corner(i), eps))
                throw std::runtime_error("Polygon point not on support plane");

        // check contains query
        if (!polygon.contains(Point(1.5*f - 1e-6*f, 0.5*f, 0.0), eps))
            throw std::runtime_error("contains() query 4 failed");
        if (!polygon.contains(Point(1.5*f, 0.5*f, 0.0), eps))
            throw std::runtime_error("contains() query 5 failed");
        if (polygon.contains(Point(1.5*f + 1e-6*f, 0.5*f, 0.0), eps))
            throw std::runtime_error("contains() query 6 failed");

        if (!polygon.contains(Point(1.25*f, 0.25*f, 0.0), eps))
            throw std::runtime_error("contains() query 7 failed");
        if (!polygon.contains(Point(1.25*f - 1e-6*f, 0.25*f, 0.0), eps))
            throw std::runtime_error("contains() query 8 failed");
        if (polygon.contains(Point(1.25*f + 1e-6*f, 0.25*f, 0.0), eps))
            throw std::runtime_error("contains() query 9 failed");

        if (!polygon.isConvex())
            throw std::runtime_error("Wrong convexity of polygon");

        // make concave polygon
        Polygon polygon2(std::vector<Point>{{Point(0.0, 0.0, 0.0),
                                             Point(f,   0.0, 0.0),
                                             Point(0.5*f, 0.5*f,   0.0),
                                             Point(f,   f,   0.0),
                                             Point(0.0, f,   0.0)}});

        if (polygon2.isConvex())
            throw std::runtime_error("Wrong convexity of polygon2");
    }

    std::cout << "All tests passed" << std::endl;
    return 0;
}
