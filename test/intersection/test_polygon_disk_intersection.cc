#include <vector>
#include <stdexcept>
#include <type_traits>

#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/polygon.hh>
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

// make quadrilateral instance
template<class Quad, class Point>
Quad makeQuad(const Point& p1, const Point& p2,
              const Point& p3, const Point& p4, std::false_type)
{ return Quad(p1, p2, p3, p4); }

// make polygon instance
template<class Quad, class Point>
Quad makeQuad(const Point& p1, const Point& p2,
              const Point& p3, const Point& p4, std::true_type)
{ return Quad(std::vector<Point>({p1, p2, p4, p3})); }

//! test quadrilateral/polygon - disk intersections
template<class QuadGeom>
void doTest()
{
    using Quad = QuadGeom;
    using ctype = typename Quad::ctype;
    using Disk = Frackit::Disk<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;
    constexpr bool isPolygonImpl = std::is_same_v<Quad, Frackit::Polygon<ctype, 3>>;

    std::vector<ctype> scales({1.0e-5, 1, 1e5});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;

        auto quad = makeQuad<Quad>(Point(-0.5*f, -0.5*f, 0.0),
                                   Point(0.5*f, -0.5*f, 0.0),
                                   Point(-0.5*f, 0.5*f, 0.0),
                                   Point(0.5*f, 0.5*f, 0.0),
                                   std::integral_constant<bool, isPolygonImpl>());

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

}

int main()
{
    using ctype = double;
    doTest<Frackit::Quadrilateral<ctype, 3>>();
    doTest<Frackit::Polygon<ctype, 3>>();
    return 0;
}
