#include <stdexcept>
#include <string>

#include <frackit/geometry/vector.hh>
#include <frackit/geometry/cylindricalsurface.hh>

//! test some functionality of circles
int main()
{
    using ctype = double;
    using CylindricalSurface = Frackit::CylindricalSurface<ctype>;
    using Point = typename CylindricalSurface::Point;
    using Vector = Frackit::Vector<ctype, 3>;

    CylindricalSurface surface(1.0, 2.0);

    // check normal vector
    using std::abs;
    if ( abs(abs(Vector(0.0, 0.0, 1.0)*Vector(surface.direction())) - 1.0) > 1e-7 )
        throw std::runtime_error(std::string("Unexpected direction"));

    // test contains() functionality
    if (surface.contains(Point(0.0 - 1e-6, 0.0, 0.0)))
        throw std::runtime_error(std::string("False positive contains() result"));

    if (surface.contains(Point(0.0, 0.0, 2.0 + 1e-6)))
        throw std::runtime_error(std::string("False positive contains() result"));

    if (surface.contains(Point(0.0 + 1e-6, 0.0, 0.0)))
        throw std::runtime_error(std::string("False positive contains() result"));

    if (surface.contains(Point(0.0, 0.0, 2.0 - 1e-6)))
        throw std::runtime_error(std::string("False positive contains() result"));

    if (surface.contains(Point(1.0 + 1e-6, 0.0, 1.0)))
        throw std::runtime_error(std::string("False positive contains() result"));

    if (surface.contains(Point(1.0 - 1e-6, 0.0, 1.0)))
        throw std::runtime_error(std::string("False positive contains() result"));

    if (!surface.cylinder().contains(Point(1.0 - 1e-6, 0.0, 1.0)))
        throw std::runtime_error(std::string("False positive contains() result on cylinder"));

    if (!surface.contains(Point(1.0, 0.0, 1.0)))
        throw std::runtime_error(std::string("False negative contains() result"));

    if (!surface.contains(Point(1.0, 0.0, 2.0)))
        throw std::runtime_error(std::string("False negative contains() result"));

    std::cout << "All tests passed" << std::endl;
    return 0;
}
