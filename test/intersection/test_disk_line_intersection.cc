#include <frackit/geometry/disk.hh>
#include <frackit/magnitude/length.hh>
#include <frackit/intersection/intersect.hh>

enum IntersectionType { point, segment, empty };

template<class G>
void checkResultGeometry(const G& geometry, IntersectionType expected)
{
    throw std::runtime_error("Unexpected intersection geometry");
}

template<int wd>
void checkResultGeometry(const Frackit::EmptyIntersection<wd>& empty, IntersectionType expected)
{
    std::cout << "Found empty intersection" << std::endl;
    if (expected != IntersectionType::empty)
        throw std::runtime_error("Got an unexpected empty intersection");
}

template<class CT, int wd>
void checkResultGeometry(const Frackit::Point<CT, wd>& p, IntersectionType expected)
{
    std::cout << "Found intersection point at " << p << std::endl;
    if (expected != IntersectionType::point)
        throw std::runtime_error("Got an unexpected point intersection");
}

template<class CT, int wd>
void checkResultGeometry(const Frackit::Segment<CT, wd>& segment, IntersectionType expected)
{
    std::cout << "Found intersection segment with corners "
              << segment.source() << " - " << segment.target()
              << " and length " << Frackit::computeLength(segment) << std::endl;
    if (expected != IntersectionType::segment)
        throw std::runtime_error("Got an unexpected segment intersection");
}

//! test disk-disk intersections
int main()
{
    using ctype = double;
    using Line = Frackit::Line<ctype, 3>;
    using Disk = Frackit::Disk<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;

    std::vector<ctype> scales({1.0e-5, 1.0, 1.0e5});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;

        Direction e1(Vector(1.0, 0.0, 0.0));
        Direction e2(Vector(0.0, 1.0, 0.0));
        Direction e3(Vector(0.0, 0.0, 1.0));
        Disk disk(Point(0.0, 0.0, 0.0), e1, e2, f, 0.5*f);

        // line that intersects in the center
        auto result = intersect(disk, Line(Point(0.0, 0.0, 0.0), e3));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);

        // line that is just outside of the disk (at major axis)
        result = intersect(disk, Line(Point(f+1e-6*f, 0.0, 0.0), e3));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::empty); }, result);

        // line that is just inside of the disk (at major axis)
        result = intersect(disk, Line(Point(f-1e-6*f, 0.0, 0.0), e3));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);

        // line that is on the rim of the disk (at major axis)
        result = intersect(disk, Line(Point(f, 0.0, 0.0), e2));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);

        // line that is just outside of the disk (at minor axis)
        result = intersect(disk, Line(Point(0.0, 0.5*f+0.5*1e-6*f, 0.0), e3));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::empty); }, result);

        // line that is just inside of the disk (at minor axis) - test flipped args
        result = intersect(Line(Point(0.0, 0.5*f-0.5*1e-6*f, 0.0), e3), disk);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);

        // line that is on the rim of the disk (at minor axis)
        result = intersect(disk, Line(Point(0.0, 0.5*f, 0.0), e1));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);

        // line that intersects in a segment on x-axis
        result = intersect(disk, Line(Point(0.0, 0.0, 0.0), e1));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result);

        std::cout << "All tests passed"  << std::endl;
    }

    return 0;
}
