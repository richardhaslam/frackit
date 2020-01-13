#include <stdexcept>
#include <string>

#include <frackit/geometry/disk.hh>

//! test some functionality of disks/ellipses
int main()
{
    using ctype = double;

    using Disk = Frackit::Disk<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;

    Vector e1(1.0, 0.0, 0.0);
    Vector e2(0.0, 1.0, 0.0);
    Disk disk(Point(0.0, 0.0, 0.0), e1, e2, 1.0, 0.5);

    // check normal vectors
    using std::abs;
    if (abs(Vector(0.0, 0.0, 1.0)*Vector(disk.normal()) - 1.0) > 1e-7)
        throw std::runtime_error(std::string("Unexpected disk normal"));

    // check point contained predicates
    if (disk.contains(Point(1.0+1e-6, 0.0, 0.0)))
        throw std::runtime_error(std::string("False positive on contained point"));

    if (!disk.contains(Point(1.0, 0.0, 0.0)))
        throw std::runtime_error(std::string("False negative on contained point"));

    if (!disk.contains(Point(0.5, 0.25, 0.0)))
        throw std::runtime_error(std::string("False negative on contained point"));

    const auto& e = disk.boundingEllipse();
    if (!e.contains(Point(1.0, 0.0, 0.0)))
        throw std::runtime_error(std::string("False negative on contained point on ellipse"));

    if (e.contains(Point(1.0-1e-6, 0.0, 0.0)))
        throw std::runtime_error(std::string("False positive on contained point on ellipse"));

    std::cout << "All tests passed" << std::endl;
    return 0;
}
