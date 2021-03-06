#include <cmath>
#include <string>
#include <stdexcept>
#include <memory>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/cylinder.hh>
#include <frackit/entitynetwork/constraints.hh>

//! test the constraints for entity networks of disks
int main()
{
    using ctype = double;

    using Disk = Frackit::Disk<ctype>;
    using Quad = Frackit::Quadrilateral<ctype, 3>;
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
        throw std::runtime_error("Did not detect distance violation");
    std::cout << "Test 1 passed" << std::endl;

    // just doesn't violates distance constraint
    Disk disk2(Point(0.0, 0.0, 0.1+1e-3), e1, e2, 1.0, 1.0);
    if (!constraints.evaluate(mainDisk, disk2))
        throw std::runtime_error("Detected false positive distance violation");
    std::cout << "Test 2 passed" << std::endl;

    // violates distance constraint (orthogonal)
    Disk disk3(Point(0.0, 0.0, 1.1-1e-3), e1, e3, 1.0, 1.0);
    if (constraints.evaluate(mainDisk, disk3))
        throw std::runtime_error("Did not detect distance violation");
    std::cout << "Test 3 passed" << std::endl;

    // just doesn't violates distance constraint (orthogonal)
    Disk disk4(Point(0.0, 0.0, 1.1+1e-3), e1, e3, 1.0, 1.0);
    if (!constraints.evaluate(mainDisk, disk4))
        throw std::runtime_error("Detected false positive distance violation");
    std::cout << "Test 4 passed" << std::endl;

    // too small intersection
    Disk disk5(Point(0.0, 0.0, 1.0-1e-4), e1, e3, 1.0, 1.0);
    if (constraints.evaluate(mainDisk, disk5))
        throw std::runtime_error("Did not detect intersection magnitude violation");
    std::cout << "Test 5 passed" << std::endl;

    // intersection magnitude ok
    Disk disk6(Point(0.0, 0.0, 0.8), e1, e3, 1.0, 1.0);
    if (!constraints.evaluate(mainDisk, disk6))
        throw std::runtime_error("False positive intersection magnitude violation");
    std::cout << "Test 6 passed" << std::endl;

    // intersection magnitude ok, but distance to boundary is too small
    Disk disk7(Point(0.38, 0.0, 0.8), e1, e3, 1.0, 1.0);
    if (constraints.evaluate(mainDisk, disk7))
        throw std::runtime_error("Did not detect intersection distance violation");
    std::cout << "Test 7 passed" << std::endl;

    // intersection angle just ok
    const Vector e22(0.0, 1.0, 1.0 + 1e-3);
    Disk disk8(Point(0.0, 0.0, 0.0), e1, e22, 1.0, 1.0);
    if (!constraints.evaluate(mainDisk, disk8))
        throw std::runtime_error("False positive intersection angle violation");
    std::cout << "Test 8 passed" << std::endl;

    // intersection angle too small
    const Vector e23(0.0, 1.0, 1.0 - 1e-3);
    Disk disk9(Point(0.0, 0.0, 0.0), e1, e23, 1.0, 1.0);
    if (constraints.evaluate(mainDisk, disk9))
        throw std::runtime_error("Did not detect intersection angle violation");
    std::cout << "Test 9 passed" << std::endl;

    // test both intersection angles also with quadrilaterals
    {
        // intersection angle just ok
        Quad quad1(Point(-1.0, -1.0, -1.0 - 0.5e-3),
                   Point( 1.0, -1.0, -1.0 - 0.5e-3),
                   Point(-1.0,  1.0, 1.0 + 0.5e-3),
                   Point( 1.0,  1.0, 1.0 + 0.5e-3));
        if (!constraints.evaluate(mainDisk, quad1))
            throw std::runtime_error("False positive intersection angle violation");
        std::cout << "Test 8 with quad passed" << std::endl;

        // intersection angle too small
        Quad quad2(Point(-1.0, -1.0, -1.0 + 0.5e-3),
                   Point( 1.0, -1.0, -1.0 + 0.5e-3),
                   Point(-1.0,  1.0, 1.0 - 0.5e-3),
                   Point( 1.0,  1.0, 1.0 - 0.5e-3));
        if (constraints.evaluate(mainDisk, quad2))
            throw std::runtime_error("Did not detect intersection angle violation");
        std::cout << "Test 9 with quad passed" << std::endl;
    }

    // Test constraints w.r.t. cylinder
    Frackit::Cylinder<ctype> cylinder(0.5, 1.0);
    Disk disk10(Point(0.0, 0.0, 0.5), e1, e2, 2.0, 2.0);
    if (!constraints.evaluate(cylinder.lateralFace(), disk10))
        throw std::runtime_error("False positive intersection distance violation");
    std::cout << "Test 10 passed" << std::endl;

    // violates intersection distance constraint
    Disk disk11(Point(0.0, 0.0, 0.951), e1, e2, 2.0, 2.0);
    if (constraints.evaluate(cylinder.lateralFace(), disk11))
        throw std::runtime_error("Did not detect intersection distance violation");
    std::cout << "Test 11 passed" << std::endl;

    // just doesn't violate intersection distance constraint
    Disk disk12(Point(0.0, 0.0, 0.849), e1, e2, 2.0, 2.0);
    if (!constraints.evaluate(cylinder.lateralFace(), disk12))
        throw std::runtime_error("False positive intersection distance violation");
    std::cout << "Test 12 passed" << std::endl;

    ///////////////////////////////////////////
    // test overloads for vectors of entities

    // 1. only admissible disks
    std::vector<Disk> admissibles({disk2, disk4, disk6, disk8});
    if (!constraints.evaluate(admissibles, mainDisk))
        throw std::runtime_error("False negative constraints violation");
    std::cout << "Test 13.1 passed" << std::endl;
    if (!constraints.evaluate(mainDisk, admissibles))
        throw std::runtime_error("False negative constraints violation");
    std::cout << "Test 13.2 passed" << std::endl;

    // 2. only non-admissible disks
    std::vector<Disk> nonadmissibles({disk1, disk3, disk5, disk7, disk9});
    if (constraints.evaluate(nonadmissibles, mainDisk))
        throw std::runtime_error("Did not detect constraints violation");
    std::cout << "Test 14.1 passed" << std::endl;
    if (constraints.evaluate(mainDisk, nonadmissibles))
        throw std::runtime_error("Did not detect constraints violation");
    std::cout << "Test 14.2 passed" << std::endl;

    // 3. mixed
    std::vector<Disk> mixed({disk1, disk2, disk3, disk4, disk5});
    if (constraints.evaluate(mixed, mainDisk))
        throw std::runtime_error("Did not detect constraints violation");
    std::cout << "Test 15.1 passed" << std::endl;
    if (constraints.evaluate(mainDisk, mixed))
        throw std::runtime_error("Did not detect constraints violation");
    std::cout << "Test 15.2 passed" << std::endl;

    ///////////////////////////////////////////
    // test overloads for pointers to abstract base class
    std::shared_ptr<Frackit::Geometry> basePtrDisk1 = std::make_shared<Disk>(disk1);
    std::shared_ptr<Frackit::Geometry> basePtrDisk2 = std::make_shared<Disk>(disk2);
    std::shared_ptr<Frackit::Geometry> basePtrMainDisk = std::make_shared<Disk>(mainDisk);

    if (constraints.evaluate(basePtrDisk1, mainDisk))
        throw std::runtime_error("Did not detect constraints violation");
    std::cout << "Test 16.1 passed" << std::endl;
    if (constraints.evaluate(mainDisk, basePtrDisk1))
        throw std::runtime_error("Did not detect constraints violation");
    std::cout << "Test 16.2 passed" << std::endl;
    if (constraints.evaluate(basePtrDisk1, basePtrMainDisk))
        throw std::runtime_error("Did not detect constraints violation");
    std::cout << "Test 16.3 passed" << std::endl;
    if (constraints.evaluate(basePtrMainDisk, basePtrDisk1))
        throw std::runtime_error("Did not detect constraints violation");
    std::cout << "Test 16.4 passed" << std::endl;

    if (!constraints.evaluate(basePtrDisk2, mainDisk))
        throw std::runtime_error("False negative constraints violation");
    std::cout << "Test 17.1 passed" << std::endl;
    if (!constraints.evaluate(mainDisk, basePtrDisk2))
        throw std::runtime_error("False negative constraints violation");
    std::cout << "Test 17.2 passed" << std::endl;
    if (!constraints.evaluate(basePtrDisk2, basePtrMainDisk))
        throw std::runtime_error("False negative constraints violation");
    std::cout << "Test 17.3 passed" << std::endl;
    if (!constraints.evaluate(basePtrMainDisk, basePtrDisk2))
        throw std::runtime_error("False negative constraints violation");
    std::cout << "Test 17.4 passed" << std::endl;

    ///////////////////////////////////////////
    // test overloads for two sets of entities

    std::vector<Disk> mainDisks({mainDisk, mainDisk});
    std::vector<std::shared_ptr<Frackit::Geometry>> mainDiskBasePtrs({basePtrMainDisk, basePtrMainDisk});

    // 1. only admissible disks
    if (!constraints.evaluate(admissibles, mainDisks))
        throw std::runtime_error("False negative constraints violation");
    std::cout << "Test 18.1 passed" << std::endl;
    if (!constraints.evaluate(mainDisks, admissibles))
        throw std::runtime_error("False negative constraints violation");
    std::cout << "Test 18.2 passed" << std::endl;
    if (!constraints.evaluate(admissibles, mainDiskBasePtrs))
        throw std::runtime_error("False negative constraints violation");
    std::cout << "Test 18.3 passed" << std::endl;
    if (!constraints.evaluate(mainDiskBasePtrs, admissibles))
        throw std::runtime_error("False negative constraints violation");
    std::cout << "Test 18.4 passed" << std::endl;

    // 2. only non-admissible disks
    if (constraints.evaluate(nonadmissibles, mainDisks))
        throw std::runtime_error("Did not detect constraints violation");
    std::cout << "Test 19.1 passed" << std::endl;
    if (constraints.evaluate(mainDisks, nonadmissibles))
        throw std::runtime_error("Did not detect constraints violation");
    std::cout << "Test 19.2 passed" << std::endl;
    if (constraints.evaluate(nonadmissibles, mainDiskBasePtrs))
        throw std::runtime_error("Did not detect constraints violation");
    std::cout << "Test 19.3 passed" << std::endl;
    if (constraints.evaluate(mainDiskBasePtrs, nonadmissibles))
        throw std::runtime_error("Did not detect constraints violation");
    std::cout << "Test 19.4 passed" << std::endl;

    // 3. mixed
    if (constraints.evaluate(mixed, mainDisks))
        throw std::runtime_error("Did not detect constraints violation");
    std::cout << "Test 20.1 passed" << std::endl;
    if (constraints.evaluate(mainDisks, mixed))
        throw std::runtime_error("Did not detect constraints violation");
    std::cout << "Test 20.2 passed" << std::endl;
    if (constraints.evaluate(mixed, mainDiskBasePtrs))
        throw std::runtime_error("Did not detect constraints violation");
    std::cout << "Test 20.3 passed" << std::endl;
    if (constraints.evaluate(mainDiskBasePtrs, mixed))
        throw std::runtime_error("Did not detect constraints violation");
    std::cout << "Test 20.4 passed" << std::endl;

    std::cout << "All tests passed" << std::endl;
    return 0;
}
