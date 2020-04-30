#include <stdexcept>
#include <string>

#include <frackit/geometry/vector.hh>
#include <frackit/geometry/hollowcylinder.hh>

//! test some functionality of circles
int main()
{
    using ctype = double;
    using HollowCylinder = Frackit::HollowCylinder<ctype>;
    using Point = HollowCylinder::Point;
    using Vector = Frackit::Vector<ctype, 3>;

    HollowCylinder hollowCylinder(0.8, 1.0, 2.0);

    // check normal vector
    using std::abs;
    if ( abs(abs(Vector(0.0, 0.0, 1.0)*Vector(hollowCylinder.axis())) - 1.0) > 1e-7 )
        throw std::runtime_error("Unexpected direction");
    std::cout << "Test 1 passed" << std::endl;

    // test contains() functionality
    if (hollowCylinder.contains(Point(0.9, 0.0, 2.0 + 1e-6)))
        throw std::runtime_error("False positive contains() result");
    std::cout << "Test 2 passed" << std::endl;

    if (!hollowCylinder.contains(Point(0.9, 0.0, 2.0 - 1e-6)))
        throw std::runtime_error("False negative contains() result");
    std::cout << "Test 3 passed" << std::endl;

    if (hollowCylinder.contains(Point(1.0 + 1e-6, 0.0, 1.0)))
        throw std::runtime_error("False positive contains() result");
    std::cout << "Test 4 passed" << std::endl;

    if (!hollowCylinder.contains(Point(1.0 - 1e-6, 0.0, 1.0)))
        throw std::runtime_error("False negative contains() result");
    std::cout << "Test 5 passed" << std::endl;

    if (hollowCylinder.contains(Point(0.8 - 1e-6, 0.0, 1.0)))
        throw std::runtime_error("False positive contains() result");
    std::cout << "Test 6 passed" << std::endl;

    if (!hollowCylinder.contains(Point(0.8 + 1e-6, 0.0, 1.0)))
        throw std::runtime_error("False negative contains() result");
    std::cout << "Test 7 passed" << std::endl;

    if (!hollowCylinder.contains(Point(1.0, 0.0, 1.0)))
        throw std::runtime_error("False negative contains() result");
    std::cout << "Test 8 passed" << std::endl;

    if (!hollowCylinder.contains(Point(1.0, 0.0, 2.0)))
        throw std::runtime_error("False negative contains() result");
    std::cout << "Test 9 passed" << std::endl;

    if (!hollowCylinder.innerLateralFace().contains(Point(0.8, 0.0, 1.0)))
        throw std::runtime_error("False negative contains() result on inner mantle");
    std::cout << "Test 10 passed" << std::endl;

    if (hollowCylinder.innerLateralFace().contains(Point(0.8 + 1e-6, 0.0, 1.0)))
        throw std::runtime_error("False positive contains() result on inner mantle");
    std::cout << "Test 11 passed" << std::endl;

    if (!hollowCylinder.innerLateralFace().contains(Point(0.8, 0.0, 2.0)))
        throw std::runtime_error("False negative contains() result on inner mantle");
    std::cout << "Test 12 passed" << std::endl;

    if (hollowCylinder.innerLateralFace().contains(Point(0.8, 0.0, 2.0 + 1e-6)))
        throw std::runtime_error("False positive contains() result on inner mantle");
    std::cout << "Test 13 passed" << std::endl;

    if (!hollowCylinder.outerLateralFace().contains(Point(1.0, 0.0, 1.0)))
        throw std::runtime_error("False negative contains() result on outer mantle");
    std::cout << "Test 14 passed" << std::endl;

    if (hollowCylinder.outerLateralFace().contains(Point(1.0 + 1e-6, 0.0, 1.0)))
        throw std::runtime_error("False positive contains() result on outer mantle");
    std::cout << "Test 15 passed" << std::endl;

    if (!hollowCylinder.outerLateralFace().contains(Point(1.0, 0.0, 2.0)))
        throw std::runtime_error("False negative contains() result on outer mantle");
    std::cout << "Test 16 passed" << std::endl;

    if (hollowCylinder.outerLateralFace().contains(Point(1.0, 0.0, 2.0 + 1e-6)))
        throw std::runtime_error("False positive contains() result on outer mantle");
    std::cout << "Test 17 passed" << std::endl;

    std::cout << "All tests passed" << std::endl;
    return 0;
}
