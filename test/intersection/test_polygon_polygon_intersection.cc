#include <vector>

#include <frackit/geometry/polygon.hh>
#include <frackit/magnitude/length.hh>
#include <frackit/intersection/intersect.hh>
#include <frackit/precision/precision.hh>

enum IntersectionType { point, segment, polygon, empty };

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
void checkResultGeometry(const Frackit::Polygon<CT, wd>& polygon, IntersectionType expected)
{
    std::cout << "Found intersection polygon with center " << polygon.center() << std::endl;
    if (expected != IntersectionType::polygon)
        throw std::runtime_error("Got an unexpected polygon intersection");
    std::cout << "Test passed" << std::endl;
}

//! test polygon-polygon intersections
int main()
{
    using ctype = double;
    using Polygon = Frackit::Polygon<ctype, 3>;
    using Point = typename Polygon::Point;
    using Direction = typename Frackit::Direction<ctype, 3>;
    using Vector = typename Direction::Vector;

    std::vector<ctype> scales({1.0e-5, 1, 1e5});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;

        Direction e11(Vector(1.0, 0.0, 0.0));
        Direction e12(Vector(0.0, 1.0, 0.0));
        using CornerVec = std::vector<Point>;
        Polygon poly1(CornerVec({Point(-0.5*f, -0.5*f, 0.0),
                                 Point(0.5*f, -0.5*f, 0.0),
                                 Point(0.5*f, 0.5*f, 0.0),
                                 Point(-0.5*f, 0.5*f, 0.0)}));

        // polygon that intersects poly1 in a segment fully in poly1
        Polygon poly2(CornerVec({Point(-0.25*f, 0.0, -0.25*f),
                                 Point(0.25*f,  0.0, -0.25*f),
                                 Point(0.25*f,  0.0,  0.25*f),
                                 Point(-0.25*f, 0.0,  0.25*f)}));
        auto result = intersect(poly1, poly2);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result);
        std::cout << "Test 1 passed" << std::endl;

        // polygon that intersects poly1 in a crooked segment
        Polygon poly3(CornerVec({Point(-0.25*f, -0.25*f, -0.25*f),
                                 Point(0.25*f,  -0.25*f, -0.25*f),
                                 Point(0.25*f,  0.25*f,  0.25*f),
                                 Point(-0.25*f, 0.25*f,  0.25*f)}));
        result = intersect(poly1, poly3);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result);
        std::cout << "Test 2 passed" << std::endl;

        // polygon that touches poly1 in origin
        Polygon poly4(CornerVec({Point(-0.25*f, -0.25*f, 0.25*f),
                                 Point(0.0*f,  0.0, 0.0*f),
                                 Point(0.25*f, 0.25*f,  0.25*f),
                                 Point(0.0*f,  0.0,  0.5*f)}));
        result = intersect(poly1, poly4);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);
        if (!std::get<Point>(result).isEqual(Point(0.0, 0.0, 0.0), Frackit::Precision<ctype>::confusion()*f))
            throw std::runtime_error("Unexpected touching point");
        std::cout << "Test 3 passed" << std::endl;

        // polygon that is slightly above poly1
        Polygon poly5(CornerVec({Point(-0.25*f, -0.25*f, 0.25*f),
                                 Point(0.0*f,  0.0, 1e-6*f),
                                 Point(0.25*f, 0.25*f,  0.25*f),
                                 Point(0.0*f,  0.0,  0.5*f)}));
        result = intersect(poly1, poly5);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::empty); }, result);
        std::cout << "Test 4 passed" << std::endl;

        std::cout << "All tests passed"  << std::endl;
    }

    return 0;
}
