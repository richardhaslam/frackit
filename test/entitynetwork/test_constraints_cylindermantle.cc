#include <cmath>
#include <string>
#include <stdexcept>

#include <frackit/geometry/disk.hh>
#include <frackit/geometry/cylindermantle.hh>
#include <frackit/entitynetwork/constraints.hh>

//! test the constraints for entity networks of disks
int main()
{
    using ctype = double;

    using CylinderMantle = Frackit::CylinderMantle<ctype>;
    using Disk = Frackit::Disk<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;

    // Basis Vectors
    const Vector e1(1.0, 0.0, 0.0);
    const Vector e2(0.0, 1.0, 0.0);
    const Vector e3(0.0, 0.0, 1.0);

    // Define constraints
    auto constraints = Frackit::makeDefaultConstraints<ctype>();
    constraints.setMinDistance(0.1);
    constraints.setMinIntersectingAngle(M_PI/4.0);
    constraints.setMinIntersectionMagnitude(0.05);
    constraints.setMinIntersectionDistance(0.05);

    // The cylinder surface to test against
    CylinderMantle cylMantle(0.5, 1.0);

    // violates distance constraint
    Disk disk1(Point(0.0, 0.0, 0.5), e1, e2, 0.4+1e-5, 0.25);
    if (constraints.evaluate(cylMantle, disk1))
        throw std::runtime_error("Did not detect distance violation");
    std::cout << "Test 1 passed" << std::endl;

    // just doesn't violate distance constraint
    Disk disk2(Point(0.0, 0.0, 0.5), e1, e2, 0.4-1e-5, 0.25);
    if (!constraints.evaluate(cylMantle, disk2))
        throw std::runtime_error("False positive distance violation");
    std::cout << "Test 2 passed" << std::endl;

    // violates the intersection magnitude constraint
    Disk disk3(Point(0.4, 0.0, 0.5), e1, e2, 0.11, 0.05);
    if (constraints.evaluate(cylMantle, disk3))
        throw std::runtime_error("Did not detect intersection magnitude violation");
    std::cout << "Test 3 passed" << std::endl;

    // does not violate the intersection magnitude constraint
    Disk disk4(Point(0.4, 0.0, 0.5), e1, e2, 0.5, 0.5);
    if (!constraints.evaluate(cylMantle, disk4))
        throw std::runtime_error("False positive intersection magnitude violation");
    std::cout << "Test 4 passed" << std::endl;

    // violates the distance to boundary constraint
    Disk disk5(Point(0.4, 0.0, 0.95 + 1e-6), e1, e2, 0.5, 0.5);
    if (constraints.evaluate(cylMantle, disk5))
        throw std::runtime_error("Did not detect distance to boundary violation");
    std::cout << "Test 5 passed" << std::endl;

    // does not violate the distance to boundary constraint
    Disk disk6(Point(0.4, 0.0, 0.95 - 1e-6), e1, e2, 0.5, 0.5);
    if (!constraints.evaluate(cylMantle, disk6))
        throw std::runtime_error("False positive intersection magnitude violation");
    std::cout << "Test 6 passed" << std::endl;

    // violates the intersection angle constraint
    const Vector e12(1.0, 0.0, 1.0 + 1e-3);
    Disk disk7(Point(0.4, 0.0, 0.5), e12, e2, 0.5, 0.5);
    if (constraints.evaluate(cylMantle, disk7))
        throw std::runtime_error("Did not detect intersection angle violation");
    std::cout << "Test 7 passed" << std::endl;

    // does not violate the distance to boundary constraint
    Disk disk8(Point(0.4, 0.0, 0.95 - 1e-6), e1, e2, 0.5, 0.5);
    if (!constraints.evaluate(cylMantle, disk6))
        throw std::runtime_error("False positive intersection magnitude violation");
    std::cout << "Test 6 passed" << std::endl;

    // TODO: Test elliptical intersection

    std::cout << "All tests passed" << std::endl;
    return 0;
}
