#include <cmath>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include <frackit/geometry/disk.hh>
#include <frackit/geometry/cylindermantle.hh>

#include <frackit/distance/distance.hh>
#include <frackit/precision/precision.hh>
#include <frackit/intersection/intersect.hh>

//! test some functionality of ellipse arcs
int main()
{
    using ctype = double;
    using Disk = Frackit::Disk<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;

    // base directions
    Direction e1(Vector(1.0, 0.0, 0.0));
    Direction e2(Vector(0.0, 1.0, 0.0));
    Direction e3(Vector(0.0, 0.0, 1.0));

    using std::abs;
    using std::sqrt;
    std::vector<ctype> scales({1e-3, 1, 1e3});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;

        // distance between two parallel disks
        auto d = Frackit::computeDistance(Disk(Point(0.0, 0.0, 0.0), e1, e2, f, 0.5*f),
                                          Disk(Point(0.0, 0.0, f), e1, e2, f, 0.5*f));
        if ( abs(d - f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error("Test 1 failed");

        d = Frackit::computeDistance(Disk(Point(0.0, 0.0, 0.0), e1, e2, f, 0.5*f),
                                     Disk(Point(f, 0.0, f), e1, e2, f, 0.5*f));
        if ( abs(d - f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error("Test 2 failed");

        d = Frackit::computeDistance(Disk(Point(0.0, 0.0, 0.0), e1, e2, f, 0.5*f),
                                     Disk(Point(3.0*f, 0.0, f), e1, e2, f, 0.5*f));
        if ( abs(d - sqrt(2.0)*f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error("Test 3 failed");

        // orthogonal disks being closest in the middle
        d = Frackit::computeDistance(Disk(Point(0.0, 0.0, 0.0), e1, e2, f, 0.5*f),
                                     Disk(Point(0.0, 0.0, 2.0*f), e1, e3, f, f));
        if ( abs(d - f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error("Test 4 failed");

        d = Frackit::computeDistance(Disk(Point(0.0, 0.0, 0.0), e1, e2, f, 0.5*f),
                                     Disk(Point(0.0, 0.0, 1.0*f + 1e-3*f), e1, e3, f, f));
        if ( abs(d - 1e-3*f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error("Test 5 failed");

        d = Frackit::computeDistance(Disk(Point(0.0, 0.0, 0.0), e1, e2, f, 0.5*f),
                                     Disk(Point(0.0, 0.0, 1.0*f), e1, e3, f, 0.5*f));
        if ( abs(d - 0.5*f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error("Test 6 failed");

        // compute the distance of an ellipse arc to a disk, which
        // results from the intersection of a disk with a cylinder surface
        Frackit::CylinderMantle<ctype> cylMantle(0.5*f, f);
        Disk disk(Point(0.25*f, 0.0, 0.5*f), e1, e2, 0.5*f, 0.25*f);
        auto is = Frackit::intersect(cylMantle, disk);
        if (is.size() != 1)
            throw std::runtime_error("Unexpected number of intersections");
        if (!std::holds_alternative<Frackit::EllipseArc<ctype, 3>>(is[0]))
            throw std::runtime_error("Unexpected intersection type");

        const auto& arc = std::get<Frackit::EllipseArc<ctype, 3>>(is[0]);
        d = Frackit::computeDistance(cylMantle.upperBoundingCircle(), arc);
        if ( abs(d - 0.5*f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error("Test 7 failed");

        // do the same with an inclined disk
        Direction e22(Vector(0.0, 1.0, 1.0));
        Disk disk2(Point(0.5*f, 0.0, 0.5*f), e1, e22, 0.5*f, 0.25*f);

        is = Frackit::intersect(cylMantle, disk2);
        if (is.size() != 1)
            throw std::runtime_error("Unexpected number of intersections");
        if (!std::holds_alternative<Frackit::EllipseArc<ctype, 3>>(is[0]))
            throw std::runtime_error("Unexpected intersection type");

        const auto& arc2 = std::get<Frackit::EllipseArc<ctype, 3>>(is[0]);
        const auto eps = Frackit::Precision<ctype>::confusion()*f;
        if (!arc2.getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f), eps))
            throw std::runtime_error("Unexpected intersection arc");

        using std::max;
        const auto zMax = max(arc2.getPoint(0.0).z(), arc2.getPoint(1.0).z());
        d = Frackit::computeDistance(cylMantle.upperBoundingCircle(), arc2);
        if ( abs(d - (f- zMax)) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error("Test 8 failed");
    }

    std::cout << "All tests passed" << std::endl;
    return 0;
}
