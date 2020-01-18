#include <cmath>
#include <string>
#include <stdexcept>

#include <frackit/geometry/disk.hh>
#include <frackit/geometry/cylinder.hh>
#include <frackit/entitynetwork/constraints.hh>

//! test the constraints for entity networks of disks
int main()
{
    using ctype = double;

    using Disk = Frackit::Disk<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;

    // Basis Vectors
    const Vector e1(1.0, 0.0, 0.0);
    const Vector e2(0.0, 1.0, 0.0);
    const Vector e3(0.0, 0.0, 1.0);

    // Define constraints
    Frackit::EntityNetworkConstraints<ctype> constraints;
    constraints.setMinDistance(0.1);
    constraints.setMinIntersectingAngle(M_PI/4.0);
    constraints.setMinIntersectionMagnitude(0.05);
    constraints.setMinIntersectionDistance(0.05);

    // disk to test the others against
    Disk mainDisk(Point(0.0, 0.0, 0.0), e1, e2, 1.0, 1.0);

    // violates distance constraint
    Disk disk1(Point(0.0, 0.0, 0.1-1e-3), e1, e2, 1.0, 1.0);
    if (constraints.evaluate(mainDisk, disk1))
        throw std::runtime_error(std::string("Did not detect distance violation"));
    std::cout << "Test 1 passed" << std::endl;

    // just doesn't violates distance constraint
    Disk disk2(Point(0.0, 0.0, 0.1+1e-3), e1, e2, 1.0, 1.0);
    if (!constraints.evaluate(mainDisk, disk2))
        throw std::runtime_error(std::string("Detected false positive distance violation"));
    std::cout << "Test 2 passed" << std::endl;

    // violates distance constraint (orthogonal)
    Disk disk3(Point(0.0, 0.0, 1.1-1e-3), e1, e3, 1.0, 1.0);
    if (constraints.evaluate(mainDisk, disk3))
        throw std::runtime_error(std::string("Did not detect distance violation"));
    std::cout << "Test 3 passed" << std::endl;

    // just doesn't violates distance constraint (orthogonal)
    Disk disk4(Point(0.0, 0.0, 1.1+1e-3), e1, e3, 1.0, 1.0);
    if (!constraints.evaluate(mainDisk, disk4))
        throw std::runtime_error(std::string("Detected false positive distance violation"));
    std::cout << "Test 4 passed" << std::endl;

    // too small intersection
    Disk disk5(Point(0.0, 0.0, 1.0-1e-4), e1, e3, 1.0, 1.0);
    if (constraints.evaluate(mainDisk, disk5))
        throw std::runtime_error(std::string("Did not detect intersection magnitude violation"));
    std::cout << "Test 5 passed" << std::endl;

    // intersection magnitude ok
    Disk disk6(Point(0.0, 0.0, 0.8), e1, e3, 1.0, 1.0);
    if (!constraints.evaluate(mainDisk, disk6))
        throw std::runtime_error(std::string("False positive intersection magnitude violation"));
    std::cout << "Test 6 passed" << std::endl;

    // intersection magnitude ok, but distance to boundary is too small
    Disk disk7(Point(0.38, 0.0, 0.8), e1, e3, 1.0, 1.0);
    if (constraints.evaluate(mainDisk, disk7))
        throw std::runtime_error(std::string("Did not detect intersection distance violation"));
    std::cout << "Test 7 passed" << std::endl;

    // intersection angle just ok
    const Vector e22(0.0, 1.0, 1.0 + 1e-3);
    Disk disk8(Point(0.0, 0.0, 0.0), e1, e22, 1.0, 1.0);
    if (!constraints.evaluate(mainDisk, disk8))
        throw std::runtime_error(std::string("False positive intersection angle violation"));
    std::cout << "Test 8 passed" << std::endl;

    // intersection angle too small
    const Vector e23(0.0, 1.0, 1.0 - 1e-3);
    Disk disk9(Point(0.0, 0.0, 0.0), e1, e23, 1.0, 1.0);
    if (constraints.evaluate(mainDisk, disk9))
        throw std::runtime_error(std::string("Did not detect intersection angle violation"));
    std::cout << "Test 9 passed" << std::endl;

    // Test constraints w.r.t. cylinder
    Frackit::Cylinder<ctype> cylinder(0.5, 1.0);
    Disk disk10(Point(0.0, 0.0, 0.5), e1, e2, 2.0, 2.0);
    if (!constraints.evaluate(cylinder.lateralFace(), disk10))
        throw std::runtime_error(std::string("False positive intersection distance violation"));
    std::cout << "Test 10 passed" << std::endl;

    // violates intersection distance constraint
    Disk disk11(Point(0.0, 0.0, 0.951), e1, e2, 2.0, 2.0);
    if (constraints.evaluate(cylinder.lateralFace(), disk11))
        throw std::runtime_error(std::string("Did not detect intersection distance violation"));
    std::cout << "Test 11 passed" << std::endl;

    // just doesn't violate intersection distance constraint
    Disk disk12(Point(0.0, 0.0, 0.849), e1, e2, 2.0, 2.0);
    if (!constraints.evaluate(cylinder.lateralFace(), disk12))
        throw std::runtime_error(std::string("False positive intersection distance violation"));
    std::cout << "Test 12 passed" << std::endl;

    std::cout << "All tests passed" << std::endl;
    return 0;
}
