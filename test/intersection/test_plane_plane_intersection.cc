#include <frackit/geometry/plane.hh>
#include <frackit/intersection/intersect.hh>

enum IntersectionType { line, plane, empty };

template<class G>
void checkResultGeometry(const G& geometry, IntersectionType expected)
{
    throw std::runtime_error(std::string("Unexpected intersection geometry"));
}

void checkResultGeometry(const Frackit::EmptyIntersection<3>& empty, IntersectionType expected)
{
    std::cout << "Found empty intersection" << std::endl;
    if (expected != IntersectionType::empty)
        throw std::runtime_error(std::string("Got an unexpected empty intersection"));
}

template<class CT>
void checkResultGeometry(const Frackit::Line<CT, 3>& line, IntersectionType expected)
{
    std::cout << "Found intersection line with support " << line.supportingPoint() << " "
              << "and direction " << Frackit::Vector<CT, 3>(line.direction()) << std::endl;
    if (expected != IntersectionType::line)
        throw std::runtime_error(std::string("Got an unexpected line intersection"));
}

template<class CT>
void checkResultGeometry(const Frackit::Plane<CT, 3>& plane, IntersectionType expected)
{
    std::cout << "Found intersection plane with support " << plane.supportingPoint() << ", "
              << "base1: " << Frackit::Vector<CT, 3>(plane.base1()) << ", "
              << "base2: " << Frackit::Vector<CT, 3>(plane.base2()) << ", "
              << "normal: " << Frackit::Vector<CT, 3>(plane.normal()) << std::endl;
    if (expected != IntersectionType::plane)
        throw std::runtime_error(std::string("Got an unexpected plane intersection"));
}

template<class CT>
void checkLineOrientation(const Frackit::Line<CT, 3>& line)
{
    const auto v1 = Frackit::Vector<CT, 3>(line.direction());
    const auto v2 = Frackit::Vector<CT, 3>(1.0, 0.0, 0.0);

    using std::abs;
    if ( abs(v1*v2) > 1e-7 )
        throw std::runtime_error(std::string("Unexpected line orientation"));
}

template<class G>
void checkLineOrientation(const G& geometry)
{
    throw std::runtime_error(std::string("Line orientation check called with non-line geometry"));
}

//! test plane-plane intersections
int main()
{
    using ctype = double;
    using Plane = Frackit::Plane<ctype, 3>;
    using Point = typename Plane::Point;
    using Direction = typename Plane::Direction;
    using Vector = typename Direction::Vector;

    std::vector<ctype> scales({1.0e-20, 1e-10, 1e-5, 1, 1e5, 1e10, 1e20});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;

        Direction e1(Vector(1.0, 0.0, 0.0));
        Direction e3(Vector(0.0, 0.0, 1.0));
        Plane plane1(Point(1.0*f, 1.0*f, 1.0*f), e3);

        // plane that intersects in the y-axis
        auto result = intersect(plane1, Plane(Point(1.0*f, 1.0*f, 1.0*f), e1));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::line); }, result);
        std::visit([&] (auto&& is) { checkLineOrientation(is); }, result);

        // plane that does not intersect
        result = intersect(plane1, Plane(Point(2.0*f, 2.0*f, 2.0*f), e3));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::empty); }, result);

        // "identical" plane
        auto sp = plane1.supportingPoint();
        auto b1 = Vector(plane1.base1());
        auto b2 = Vector(plane1.base2());
        b1 *= f; b2 *= f;
        sp += b1; sp += b2;

        result = intersect(plane1, Plane(sp, e3));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::plane); }, result);

        std::cout << "All tests passed"  << std::endl;
    }

    return 0;
}
