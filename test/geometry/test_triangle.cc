#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

#include <frackit/geometry/triangle.hh>
#include <frackit/precision/precision.hh>

//! test some functionality of triangles
int main()
{
    using ctype = double;

    using Triangle = Frackit::Triangle<ctype, 3>;
    using Point = typename Triangle::Point;

    std::vector<ctype> scales({1e-5, 1, 1e5});
    for (auto f : scales)
    {
        Triangle triangle(Point(0.0, 0.0, 0.0),
                          Point(f,   0.0, 0.0),
                          Point(0.0, f,   0.0));
        using std::abs;
        using std::sqrt;
        const auto eps = Frackit::Precision<ctype>::confusion()*f;
        if ( abs(triangle.edge(0).length() - f) > eps )
            throw std::runtime_error(std::string("Unexpected length of edge 0"));
        if ( abs(triangle.edge(1).length() - sqrt(2.0)*f) > eps )
            throw std::runtime_error(std::string("Unexpected length of edge 1"));
        if ( abs(triangle.edge(2).length() - f) > eps )
            throw std::runtime_error(std::string("Unexpected length of edge 2"));
        if ( abs(triangle.area() - 0.5*f*f) > eps*eps )
            throw std::runtime_error(std::string("Unexpected triangle area"));

        // check if points are on supporting plane
        const auto plane = triangle.supportingPlane();
        if (!plane.contains(triangle.corner(0), eps))
            throw std::runtime_error(std::string("Point 0 not on support plane"));
        if (!plane.contains(triangle.corner(1), eps))
            throw std::runtime_error(std::string("Point 1 not on support plane"));
        if (!plane.contains(triangle.corner(2), eps))
            throw std::runtime_error(std::string("Point 2 not on support plane"));
        if (plane.contains(Point(0.0, 0.0, 1e-6*f), eps))
            throw std::runtime_error(std::string("Point should not be on plane"));

        // check contains() query
        if (!triangle.contains(Point(0.5*f, 0.0, 0.0), eps))
            throw std::runtime_error(std::string("contains() query 1 failed"));
        if (!triangle.contains(Point(0.5*f, 1e-6*f, 0.0), eps))
            throw std::runtime_error(std::string("contains() query 2 failed"));
        if (triangle.contains(Point(0.5*f, -1e-6*f, 0.0), eps))
            throw std::runtime_error(std::string("contains() query 3 failed"));
    }

    std::cout << "All tests passed" << std::endl;
    return 0;
}
