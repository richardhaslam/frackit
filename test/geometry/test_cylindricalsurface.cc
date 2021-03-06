#include <stdexcept>
#include <string>

#include <frackit/geometry/vector.hh>
#include <frackit/geometry/cylindersurface.hh>

//! test some functionality of circles
int main()
{
    using ctype = double;
    using CylinderSurface = Frackit::CylinderSurface<ctype>;
    using Point = typename CylinderSurface::Point;
    using Vector = Frackit::Vector<ctype, 3>;

    CylinderSurface surface(1.0, 2.0);

    // check normal vector
    using std::abs;
    if ( abs(abs(Vector(0.0, 0.0, 1.0)*Vector(surface.direction())) - 1.0) > 1e-7 )
        throw std::runtime_error("Unexpected direction");

    // test contains() functionality
    if (surface.contains(Point(0.0 - 1e-6, 0.0, 0.0)))
        throw std::runtime_error("False positive contains() result");

    if (surface.contains(Point(0.0, 0.0, 2.0 + 1e-6)))
        throw std::runtime_error("False positive contains() result");

    if (surface.contains(Point(0.0 + 1e-6, 0.0, 0.0)))
        throw std::runtime_error("False positive contains() result");

    if (surface.contains(Point(0.0, 0.0, 2.0 - 1e-6)))
        throw std::runtime_error("False positive contains() result");

    if (surface.contains(Point(1.0 + 1e-6, 0.0, 1.0)))
        throw std::runtime_error("False positive contains() result");

    if (surface.contains(Point(1.0 - 1e-6, 0.0, 1.0)))
        throw std::runtime_error("False positive contains() result");

    if (!surface.cylinder().contains(Point(1.0 - 1e-6, 0.0, 1.0)))
        throw std::runtime_error("False positive contains() result on cylinder");

    if (!surface.contains(Point(1.0, 0.0, 1.0)))
        throw std::runtime_error("False negative contains() result");

    if (!surface.contains(Point(1.0, 0.0, 2.0)))
        throw std::runtime_error("False negative contains() result");

    std::cout << "All tests passed" << std::endl;
    return 0;
}
