#include <vector>
#include <TopExp.hxx>

#include <frackit/geometry/disk.hh>
#include <frackit/geometry/polygon.hh>
#include <frackit/magnitude/length.hh>
#include <frackit/magnitude/area.hh>
#include <frackit/precision/precision.hh>
#include <frackit/occ/breputilities.hh>

#include <frackit/intersection/intersect.hh>

enum IntersectionType { point, edge, face, empty };

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

void checkResultGeometry(const TopoDS_Edge& edge, IntersectionType expected)
{
    std::cout << "Found intersection edge with corners "
              << Frackit::OCCUtilities::point(TopExp::FirstVertex(edge)) << " - "
              << Frackit::OCCUtilities::point(TopExp::LastVertex(edge))
              << " and length " << Frackit::computeLength(edge) << std::endl;
    if (expected != IntersectionType::edge)
        throw std::runtime_error("Got an unexpected segment intersection");
}

void checkResultGeometry(const TopoDS_Face& face, IntersectionType expected)
{
    std::cout << "Found intersection face with area "
              << Frackit::computeArea(face) << std::endl;
    if (expected != IntersectionType::face)
        throw std::runtime_error("Got an unexpected face intersection");
}

//! test polygon-face intersections
int main()
{
    using ctype = double;
    using Polygon = Frackit::Polygon<ctype, 3>;
    using Point = typename Polygon::Point;
    using PointVector = std::vector<Point>;

    std::vector<ctype> scales({1.0e-3, 1.0, 1.0e3});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;
        const auto eps = Frackit::Precision<ctype>::confusion()*f;

        // Make polygon
        Polygon polyGeom(PointVector({Point(0.0, 0.0, 0.0), Point(f, 0.0, 0.0),
                                      Point(f, f,   0.0), Point(0.0, f,   0.0)}));
        const auto polyFace = Frackit::OCCUtilities::getShape(polyGeom);

        // create a polygon that intersects in an edge of length f
        Polygon poly1(PointVector({Point(-f,  0.0, -f), Point(f, 0.0, -f),
                                   Point(f, 0.0,  f), Point(-f,  0.0,  f)}));
        auto result = intersect(poly1, polyFace);
        if (result.size() != 1)
            throw std::runtime_error("Unexpected intersection size");
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::edge); }, result[0]);

        using std::abs;
        if ( abs(Frackit::computeLength(std::get<TopoDS_Edge>(result[0])) - f) > eps )
            throw std::runtime_error("Unexpected intersection edge length");
        std::cout << "Test 1 passed" << std::endl;

        // create a polygon that intersects in an edge of length 0.5*f
        Polygon poly2(PointVector({Point(-0.5*f,  0.0, -f), Point(0.5*f, 0.0, -f),
                                   Point(0.5*f, 0.0,  f), Point(-0.5*f,  0.0,  f)}));
        result = intersect(poly2, polyFace);
        if (result.size() != 1)
            throw std::runtime_error("Unexpected intersection size");
        for (const auto& isection : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::edge); }, isection);

        using std::abs;
        if ( abs(Frackit::computeLength(std::get<TopoDS_Edge>(result[0])) - 0.5*f) > eps )
            throw std::runtime_error("Unexpected intersection edge 1 length");
        std::cout << "Test 2 passed" << std::endl;

        // 1 touching point
        Polygon poly3(PointVector({Point(-0.5*f,  0.0, 0.5*f), Point(0.0,   0.0, 0.0),
                                   Point(0.5*f, 0.0, 0.5*f), Point( 0.0,    0.0,     f)}));
        result = intersect(poly3, polyFace);
        if (result.size() != 1)
            throw std::runtime_error("Unexpected intersection size");
        for (const auto& isection : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, isection);
        std::cout << "Test 3 passed" << std::endl;

        // 1 touching point (in-plane)
        Polygon poly4(PointVector({Point(-0.5*f,  0.0, 0.0), Point(0.0, 0.5*f, 0.0),
                                   Point(-0.25*f, 1.0*f, 0.0), Point(-0.5*f,  0.5*f, 0.0)}));
        result = intersect(poly4, polyFace);
        if (result.size() != 1)
            throw std::runtime_error("Unexpected intersection size");
        for (const auto& isection : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, isection);
        std::cout << "Test 4 passed" << std::endl;

        // no touching point
        Polygon poly5(PointVector({Point(-0.5*f,  0.0, 0.5*f), Point(0.0,   0.0, 1e-3*f),
                                   Point(0.5*f, 0.0, 0.5*f), Point( 0.0,    0.0,     f)}));
        result = intersect(poly5, polyFace);
        if (result.size() != 1)
            throw std::runtime_error("Unexpected intersection size");
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::empty); }, result[0]);
        std::cout << "Test 5 passed" << std::endl;

        // intersection face that has an area of f^2
        Polygon poly6(PointVector({Point(-2.0*f, -2.0*f, 0.0), Point(2.0*f, -2.0*f, 0.0),
                                   Point(2.0*f,  2.0*f, 0.0), Point(-2.0*f,  2.0*f, 0.0)}));
        result = intersect(poly6, polyFace);
        if (result.size() != 1)
            throw std::runtime_error("Unexpected intersection size");
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::face); }, result[0]);
        if ( abs(Frackit::computeArea(std::get<TopoDS_Face>(result[0])) - f*f) > eps*eps )
            throw std::runtime_error("Unexpected intersection face area");
        std::cout << "Test 5 passed" << std::endl;

        std::cout << "All tests passed"  << std::endl;
    }

    return 0;
}
