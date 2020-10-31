#include <frackit/geometry/quadrilateral.hh>
#include <frackit/magnitude/length.hh>
#include <frackit/intersection/intersect.hh>
#include <frackit/precision/precision.hh>

enum IntersectionType { point, segment, quad, empty };

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
void checkResultGeometry(const Frackit::Quadrilateral<CT, wd>& quad, IntersectionType expected)
{
    std::cout << "Found intersection quad with center " << quad.center() << std::endl;
    if (expected != IntersectionType::quad)
        throw std::runtime_error("Got an unexpected quadrilateral intersection");
    std::cout << "Test passed" << std::endl;
}

//! test quadrilateral-quadrilateral intersections
int main()
{
    using ctype = double;
    using Quad = Frackit::Quadrilateral<ctype, 3>;
    using Point = typename Quad::Point;
    using Direction = typename Frackit::Direction<ctype, 3>;
    using Vector = typename Direction::Vector;

    std::vector<ctype> scales({1.0e-5, 1, 1e5});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;

        Direction e11(Vector(1.0, 0.0, 0.0));
        Direction e12(Vector(0.0, 1.0, 0.0));
        Quad quad1(Point(-0.5*f, -0.5*f, 0.0),
                   Point(0.5*f, -0.5*f, 0.0),
                   Point(-0.5*f, 0.5*f, 0.0),
                   Point(0.5*f, 0.5*f, 0.0));

        // quad that intersects quad1 in a segment fully in quad1
        Quad quad2(Point(-0.25*f, 0.0, -0.25*f),
                   Point(0.25*f,  0.0, -0.25*f),
                   Point(-0.25*f, 0.0,  0.25*f),
                   Point(0.25*f,  0.0,  0.25*f));
        auto result = intersect(quad1, quad2);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result);
        std::cout << "Test 1 passed" << std::endl;

        // quad that intersects quad1 in a crooked segment
        Quad quad3(Point(-0.25*f, -0.25*f, -0.25*f),
                   Point(0.25*f,  -0.25*f, -0.25*f),
                   Point(-0.25*f, 0.25*f,  0.25*f),
                   Point(0.25*f,  0.25*f,  0.25*f));
        result = intersect(quad1, quad3);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result);
        std::cout << "Test 2 passed" << std::endl;

        // quad that touches quad1 in origin
        Quad quad4(Point(-0.25*f, -0.25*f, 0.25*f),
                   Point(0.0*f,  0.0, 0.0*f),
                   Point(0.0*f,  0.0,  0.5*f),
                   Point(0.25*f, 0.25*f,  0.25*f));
        result = intersect(quad1, quad4);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);
        if (!std::get<Point>(result).isEqual(Point(0.0, 0.0, 0.0), Frackit::Precision<ctype>::confusion()*f))
            throw std::runtime_error("Unexpected touching point");
        std::cout << "Test 3 passed" << std::endl;

        // quad that is slightly above quad1
        Quad quad5(Point(-0.25*f, -0.25, 0.25*f),
                   Point(0.0*f,  0.0, 1e-6*f),
                   Point(0.0*f,  0.0,  0.5*f),
                   Point(0.25*f, 0.25,  0.25*f));
        result = intersect(quad1, quad5);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::empty); }, result);
        std::cout << "Test 4 passed" << std::endl;

        std::cout << "All tests passed"  << std::endl;
    }

    return 0;
}
