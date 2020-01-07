#include <cmath>
#include <string>
#include <stdexcept>

#include <frackit/geometry/disk.hh>
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
    constraints.setMinIntersectionMagnitude(0.01);
    constraints.setMinIntersectionDistance(0.05);

    // disk to test the others against
    Disk mainDisk(Point(0.0, 0.0, 0.0), e1, e2, 1.0, 1.0);

    // violates distance constraint
    Disk disk1(Point(0.0, 0.0, 0.1-1e-6), e1, e2, 1.0, 1.0);
    if (constraints.evaluate(mainDisk, disk1))
        throw std::runtime_error(std::string("Did not detect distance violation"));

    std::cout << "All tests passed" << std::endl;
    return 0;
}
