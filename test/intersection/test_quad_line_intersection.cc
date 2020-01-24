#include <frackit/geometry/quadrilateral.hh>
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

//! test quadrilateral-quadrilateral intersections
int main()
{
    using ctype = double;
    using Line = Frackit::Line<ctype, 3>;
    using Quad = Frackit::Quadrilateral<ctype, 3>;
    using Point = typename Quad::Point;
    using Direction = typename Line::Direction;
    using Vector = typename Direction::Vector;

    std::vector<ctype> scales({1.0e-5, 1.0, 1.0e5});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;

        Direction e1(Vector(1.0, 0.0, 0.0));
        Direction e12(Vector(1.0, -1.0, 0.0));
        Direction e2(Vector(0.0, 1.0, 0.0));
        Direction e3(Vector(0.0, 0.0, 1.0));
        Quad quad(Point(-1.0*f, -0.5*f, 0.0),
                  Point(1.0*f, -0.5*f, 0.0),
                  Point(-1.0*f, 0.5*f, 0.0),
                  Point(1.0*f, 0.5*f, 0.0));

        // line that intersects in the center
        auto result = intersect(quad, Line(Point(0.0, 0.0, 0.0), e3));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);
        std::cout << "Test 1 passed" << std::endl;

        // line that is just outside of the quad (on the right side)
        result = intersect(quad, Line(Point(f+1e-6*f, 0.0, 0.0), e3));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::empty); }, result);
        std::cout << "Test 2 passed" << std::endl;

        // line that is just inside of the quad (on the right side)
        result = intersect(quad, Line(Point(f-1e-6*f, 0.0, 0.0), e3));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);
        std::cout << "Test 3 passed" << std::endl;

        // line that intersects on the rim of the quad
        result = intersect(quad, Line(Point(f, 0.0, 0.0), e3));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);
        std::cout << "Test 4 passed" << std::endl;

        // line that is just outside of the quad (top side)
        result = intersect(quad, Line(Point(0.0, 0.5*f+0.5*1e-6*f, 0.0), e3));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::empty); }, result);
        std::cout << "Test 5 passed" << std::endl;

        // line that is just inside of the quad (top side) - test flipped args
        result = intersect(Line(Point(0.0, 0.5*f-0.5*1e-6*f, 0.0), e3), quad);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);
        std::cout << "Test 6 passed" << std::endl;

        // line that is on the rim of the quad (at corner axis)
        result = intersect(quad, Line(Point(1.0*f, 0.5*f, 0.0), e12));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);
        std::cout << "Test 7 passed" << std::endl;

        // line that intersects in a segment on x-axis
        result = intersect(quad, Line(Point(0.0, 0.0, 0.0), e1));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result);
        std::cout << "Test 8 passed" << std::endl;

        // line that describes the upper part of the quad rim
        result = intersect(quad, Line(Point(0.0, 0.5*f, 0.0), e1));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result);
        std::cout << "Test 9 passed" << std::endl;

        std::cout << "All tests passed"  << std::endl;
    }

    return 0;
}
