#include <frackit/geometry/disk.hh>
#include <frackit/intersection/intersect.hh>
#include <frackit/intersection/intersectionresult.hh>

enum IntersectionType { point, segment, disk, empty };

template<class G>
void checkResultGeometry(const G& geometry, IntersectionType expected)
{
    throw std::runtime_error(std::string("Unexpected intersection geometry"));
}

void checkResultGeometry(const Frackit::EmptyIntersection<3>& geometry, IntersectionType expected)
{
    std::cout << "Found empty intersection" << std::endl;
    if (expected != IntersectionType::empty)
        throw std::runtime_error(std::string("Unexpected empty intersection"));
    std::cout << "Test passed" << std::endl;
}

template<class CT, int wd>
void checkResultGeometry(const Frackit::Point<CT, wd>& p, IntersectionType expected)
{
    std::cout << "Found intersection point at " << p << std::endl;
    if (expected != IntersectionType::point)
        throw std::runtime_error(std::string("Got an unexpected point intersection"));
    std::cout << "Test passed" << std::endl;
}

template<class CT, int wd>
void checkResultGeometry(const Frackit::Segment<CT, wd>& segment, IntersectionType expected)
{
    std::cout << "Found intersection segment with corners "
              << segment.source() << " - " << segment.target()
              << " and length " << segment.length() << std::endl;
    if (expected != IntersectionType::segment)
        throw std::runtime_error(std::string("Got an unexpected segment intersection"));
    std::cout << "Test passed" << std::endl;
}

template<class CT, int wd>
void checkResultGeometry(const Frackit::Disk<CT>& disk, IntersectionType expected)
{
    std::cout << "Found intersection disk with center " << disk.center() << ", "
              << "major axis: " << Frackit::Vector<CT, wd>(disk.majorAxis()) << ", "
              << "minor axis: " << Frackit::Vector<CT, wd>(disk.minorAxis()) << ", "
              << "major ax length: " << disk.majorAxisLength() << ", "
              << "minor ax length: " << disk.minorAxisLength() << std::endl;
    if (expected != IntersectionType::disk)
        throw std::runtime_error(std::string("Got an unexpected segment intersection"));
    std::cout << "Test passed" << std::endl;
}

//! test disk-disk intersections
int main()
{
    using ctype = double;
    using Disk = Frackit::Disk<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;

    std::vector<ctype> scales({1.0e-5, 1, 1e5});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;

        Direction e11(Vector(1.0, 0.0, 0.0));
        Direction e12(Vector(0.0, 1.0, 0.0));
        Disk disk1(Point(0.0, 0.0, 0.0), e11, e12, f, 0.5*f);

        // disk that intersects disk 1 in a segment fully in disk1
        Direction e22(Vector(0.0, 0.0, 1.0));
        auto result = intersect(disk1, Disk(Point(0.0, 0.0, 0.0), e11, e22, 0.75*f, 0.5*f));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result);

        // disk that intersects disk 1 in a crooked segment
        Direction e23(Vector(1.0, 1.0, 0.0));
        Direction e24(Vector(-1.0, 1.0, 1.0));
        result = intersect(disk1, Disk(Point(0.0, 0.0, 0.0), e23, e24, 0.75*f, 0.5*f));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result);

        // disk that touches disk 1 in origin
        result = intersect(disk1, Disk(Point(0.0, 0.0, f), e11, e22, f, f));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);

        // disk that is slightly above disk 1
        result = intersect(disk1, Disk(Point(0.0, 0.0, f+1e-6*f), e11, e22, f, f));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::empty); }, result);

        std::cout << "All tests passed"  << std::endl;
    }

    return 0;
}
