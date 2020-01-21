#include <stdexcept>
#include <string>

#include <frackit/geometry/vector.hh>
#include <frackit/geometry/cylinder.hh>

//! test some functionality of circles
int main()
{
    using ctype = double;
    using Cylinder = Frackit::Cylinder<ctype>;
    using Point = typename Cylinder::Point;
    using Vector = Frackit::Vector<ctype, 3>;

    Cylinder cylinder(1.0, 2.0);

    // check normal vector
    using std::abs;
    if ( abs(abs(Vector(0.0, 0.0, 1.0)*Vector(cylinder.direction())) - 1.0) > 1e-7 )
        throw std::runtime_error("Unexpected direction");
    std::cout << "Test 1 passed" << std::endl;

    // test contains() functionality
    if (cylinder.contains(Point(0.0, 0.0, 2.0 + 1e-6)))
        throw std::runtime_error("False positive contains() result");
    std::cout << "Test 2 passed" << std::endl;

    if (!cylinder.contains(Point(0.0, 0.0, 2.0 - 1e-6)))
        throw std::runtime_error("False negative contains() result");
    std::cout << "Test 3 passed" << std::endl;

    if (cylinder.contains(Point(1.0 + 1e-6, 0.0, 1.0)))
        throw std::runtime_error("False positive contains() result");
    std::cout << "Test 4 passed" << std::endl;

    if (!cylinder.contains(Point(1.0 - 1e-6, 0.0, 1.0)))
        throw std::runtime_error("False negative contains() result");
    std::cout << "Test 5 passed" << std::endl;

    if (!cylinder.contains(Point(1.0, 0.0, 1.0)))
        throw std::runtime_error("False negative contains() result");
    std::cout << "Test 6 passed" << std::endl;

    if (!cylinder.contains(Point(1.0, 0.0, 2.0)))
        throw std::runtime_error("False negative contains() result");

    std::cout << "All tests passed" << std::endl;
    return 0;
}
