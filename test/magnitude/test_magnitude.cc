#include <iostream>
#include <vector>
#include <stdexcept>
#include <string>

#include <frackit/magnitude/containedmagnitude.hh>
#include <frackit/precision/defaultepsilon.hh>
#include <frackit/occ/breputilities.hh>

#include <frackit/geometry/box.hh>
#include <frackit/geometry/disk.hh>

//! test computation of lengths/areas/volumes
int main()
{
    using namespace Frackit;
    using ctype = double;
    using Disk = Disk<ctype>;
    using Box = Box<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;

    // basis directions
    const Direction e1(Vector(1.0, 0.0, 0.0));
    const Direction e2(Vector(0.0, 1.0, 0.0));
    const Direction e3(Vector(0.0, 0.0, 1.0));

    std::vector<ctype> scales({1});
    for (auto f : scales)
    {
        Box box(0.0, 0.0, 0.0, f, f, f);

        // disk outside the box. Area should be zero.
        Disk disk(Point(1.5*f, 0.0, 0.5*f), e1, e2, 0.5*f, 0.5*f);
        if (computeContainedMagnitude(disk, box))
            throw std::runtime_error(std::string("Test 1 failed"));
        std::cout << "Test 1 passed" << std::endl;

        // disk area should be half of the original one.
        using std::abs;
        Disk disk2(Point(1.0*f, 0.5*f, 0.5*f), e1, e2, 0.5*f, 0.5*f);
        const auto eps = defaultEpsilon(disk2);
        if ( abs(computeContainedMagnitude(disk2, box) - 0.5*disk.area()) > eps*eps)
            throw std::runtime_error(std::string("Test 2 failed"));
        std::cout << "Test 2 passed" << std::endl;
    }

    std::cout << "All tests passed" << std::endl;

    return 0;
}
