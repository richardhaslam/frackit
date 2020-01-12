#include <TopExp.hxx>

#include <frackit/geometry/disk.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/magnitude/length.hh>
#include <frackit/magnitude/area.hh>
#include <frackit/precision/precision.hh>
#include <frackit/occ/breputilities.hh>

#include <frackit/intersection/intersect.hh>
#include <frackit/intersection/intersectionresult.hh>

enum IntersectionType { point, edge, face, empty };

template<class G>
void checkResultGeometry(const G& geometry, IntersectionType expected)
{
    throw std::runtime_error(std::string("Unexpected intersection geometry"));
}

template<int wd>
void checkResultGeometry(const Frackit::EmptyIntersection<wd>& empty, IntersectionType expected)
{
    std::cout << "Found empty intersection" << std::endl;
    if (expected != IntersectionType::empty)
        throw std::runtime_error(std::string("Got an unexpected empty intersection"));
}

template<class CT, int wd>
void checkResultGeometry(const Frackit::Point<CT, wd>& p, IntersectionType expected)
{
    std::cout << "Found intersection point at " << p << std::endl;
    if (expected != IntersectionType::point)
        throw std::runtime_error(std::string("Got an unexpected point intersection"));
}

void checkResultGeometry(const TopoDS_Edge& edge, IntersectionType expected)
{
    std::cout << "Found intersection edge with corners "
              << Frackit::OCCUtilities::point(TopExp::FirstVertex(edge)) << " - "
              << Frackit::OCCUtilities::point(TopExp::LastVertex(edge))
              << " and length " << Frackit::computeLength(edge) << std::endl;
    if (expected != IntersectionType::edge)
        throw std::runtime_error(std::string("Got an unexpected segment intersection"));
}

void checkResultGeometry(const TopoDS_Face& face, IntersectionType expected)
{
    std::cout << "Found intersection face with area "
              << Frackit::computeArea(face) << std::endl;
    if (expected != IntersectionType::face)
        throw std::runtime_error(std::string("Got an unexpected face intersection"));
}
//! test disk-face intersections
int main()
{
    using ctype = double;
    using Quad = Frackit::Quadrilateral<ctype, 3>;
    using Disk = Frackit::Disk<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;

    std::vector<ctype> scales({1.0e-3, 1.0, 1.0e3});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;
        const auto eps = Frackit::Precision<ctype>::confusion()*f;

        // Make quad
        Quad quadGeom = Quad(Point(0.0, 0.0, 0.0), Point(f, 0.0, 0.0),
                             Point(0.0, f,   0.0), Point(f, f,   0.0));
        const auto quadFace = Frackit::OCCUtilities::getShape(quadGeom);

        // basis vectors for the disks
        Direction e1(Vector(1.0, 0.0, 0.0));
        Direction e2(Vector(0.0, 1.0, 0.0));
        Direction e3(Vector(0.0, 0.0, 1.0));

        // create a disk that intersects in an edge of length 2*radius
        Disk disk(Point(0.5*f, 0.5*f, 0.0), e1, e3, 0.25*f, 0.25*f);
        auto result = intersect(disk, quadFace);
        if (result.size() != 1)
            throw std::runtime_error(std::string("Unexpected intersection size"));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::edge); }, result[0]);

        using std::abs;
        if ( abs(Frackit::computeLength(std::get<TopoDS_Edge>(result[0])) - 0.5*f) > eps )
            throw std::runtime_error(std::string("Unexpected intersection edge length"));
        std::cout << "Test 1 passed" << std::endl;

        // create a disk that intersects in an edge of length radius
        Disk disk2(Point(f, 0.5*f, 0.0), e1, e3, 0.25*f, 0.25*f);
        result = intersect(disk2, quadFace);
        if (result.size() != 1)
            throw std::runtime_error(std::string("Unexpected intersection size"));
        for (const auto& isection : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::edge); }, isection);

        using std::abs;
        if ( abs(Frackit::computeLength(std::get<TopoDS_Edge>(result[0])) - 0.25*f) > eps )
            throw std::runtime_error(std::string("Unexpected intersection edge 1 length"));
        std::cout << "Test 2 passed" << std::endl;

        // 1 touching point
        Disk disk3(Point(0.5*f, 0.5*f, 0.5*f), e1, e3, 0.5*f, 0.5*f);
        result = intersect(disk3, quadFace);
        if (result.size() != 1)
            throw std::runtime_error(std::string("Unexpected intersection size"));
        for (const auto& isection : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, isection);
        std::cout << "Test 3 passed" << std::endl;

        // 1 touching point
        Disk disk4(Point(1.5*f, 0.5*f, 0.0), e1, e2, 0.5*f, 0.5*f);
        result = intersect(disk4, quadFace);
        if (result.size() != 1)
            throw std::runtime_error(std::string("Unexpected intersection size"));
        for (const auto& isection : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, isection);
        std::cout << "Test 4 passed" << std::endl;

        // no touching point
        Disk disk5(Point(1.5*f + 1e-3*f, 0.5*f, 0.0), e1, e2, 0.5*f, 0.5*f);
        result = intersect(disk5, quadFace);
        if (result.size() != 1)
            throw std::runtime_error(std::string("Unexpected intersection size"));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::empty); }, result[0]);
        std::cout << "Test 5 passed" << std::endl;

        // intersection face that has an area of f^2
        Disk disk6(Point(0.5*f, 0.5*f, 0.0), e1, e2, 2.0*f, 2.0*f);
        result = intersect(disk6, quadFace);
        if (result.size() != 1)
            throw std::runtime_error(std::string("Unexpected intersection size"));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::face); }, result[0]);
        if ( abs(Frackit::computeArea(std::get<TopoDS_Face>(result[0])) - f*f) > eps*eps )
            throw std::runtime_error(std::string("Unexpected intersection face area"));
        std::cout << "Test 5 passed" << std::endl;

        std::cout << "All tests passed"  << std::endl;
    }

    return 0;
}
