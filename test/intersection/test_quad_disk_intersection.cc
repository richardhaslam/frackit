#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/magnitude/length.hh>
#include <frackit/intersection/intersect.hh>

enum IntersectionType { point, segment, face, empty };

template<class G>
void checkResultGeometry(const G& geometry, IntersectionType expected)
{
    throw std::runtime_error("Unexpected intersection geometry");
}

void checkResultGeometry(const Frackit::EmptyIntersection<3>& geometry, IntersectionType expected)
{
    std::cout << "Found empty intersection" << std::endl;
    if (expected != IntersectionType::empty)
        throw std::runtime_error("Unexpected empty intersection");
    std::cout << "Test passed" << std::endl;
}

template<class CT, int wd>
void checkResultGeometry(const Frackit::Point<CT, wd>& p, IntersectionType expected)
{
    std::cout << "Found intersection point at " << p << std::endl;
    if (expected != IntersectionType::point)
        throw std::runtime_error("Got an unexpected point intersection");
    std::cout << "Test passed" << std::endl;
}

template<class CT, int wd>
void checkResultGeometry(const Frackit::Segment<CT, wd>& segment, IntersectionType expected)
{
    std::cout << "Found intersection segment with corners "
              << segment.source() << " - " << segment.target()
              << " and length " << computeLength(segment) << std::endl;
    if (expected != IntersectionType::segment)
        throw std::runtime_error("Got an unexpected segment intersection");
    std::cout << "Test passed" << std::endl;
}

template<class CT, int wd>
void checkResultGeometry(const TopoDS_Face& face, IntersectionType expected)
{
    std::cout << "Found intersection face" << std::endl;
    if (expected != IntersectionType::face)
        throw std::runtime_error("Got an unexpected face intersection");
    std::cout << "Test passed" << std::endl;
}

//! test quadrilateral-disk intersections
int main()
{
    using ctype = double;
    using Disk = Frackit::Disk<ctype>;
    using Quad = Frackit::Quadrilateral<ctype, 3>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;

    std::vector<ctype> scales({1.0e-5, 1, 1e5});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;

        Quad quad(Point(-0.5*f, -0.5*f, 0.0),
                  Point(0.5*f, -0.5*f, 0.0),
                  Point(-0.5*f, 0.5*f, 0.0),
                  Point(0.5*f, 0.5*f, 0.0));

        // disk that intersects quad in a segment fully in quad
        Direction e11(Vector(1.0, 0.0, 0.0));
        Direction e22(Vector(0.0, 0.0, 1.0));
        auto result = intersect(quad, Disk(Point(0.0, 0.0, 0.0), e11, e22, 0.75*f, 0.5*f));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result);
        std::cout << "Test 1 passed" << std::endl;

        // disk that intersects quad in a crooked segment
        Direction e23(Vector(1.0, 1.0, 0.0));
        Direction e24(Vector(-1.0, 1.0, 1.0));
        result = intersect(quad, Disk(Point(0.0, 0.0, 0.0), e23, e24, 0.75*f, 0.5*f));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result);
        std::cout << "Test 2 passed" << std::endl;

        // disk that touches quad in origin
        result = intersect(quad, Disk(Point(0.0, 0.0, f), e11, e22, f, f));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);
        std::cout << "Test 3 passed" << std::endl;

        // disk that is slightly above quad
        result = intersect(quad, Disk(Point(0.0, 0.0, f+1e-6*f), e11, e22, f, f));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::empty); }, result);
        std::cout << "Test 4 passed" << std::endl;

        std::cout << "All tests passed"  << std::endl;
    }

    return 0;
}
