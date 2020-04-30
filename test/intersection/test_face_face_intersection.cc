#include <cmath>
#include <stdexcept>
#include <TopExp.hxx>

#include <frackit/geometry/disk.hh>
#include <frackit/geometry/quadrilateral.hh>
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
              << " and length " << std::setprecision(20) << Frackit::computeLength(edge) << std::endl;
    if (expected != IntersectionType::edge)
        throw std::runtime_error("Got an unexpected segment intersection");
}

void checkResultGeometry(const TopoDS_Face& face, IntersectionType expected)
{
    std::cout << "Found intersection face with area "
              << std::setprecision(20) << Frackit::computeArea(face) << std::endl;
    if (expected != IntersectionType::face)
        throw std::runtime_error("Got an unexpected face intersection");
}
//! test disk-face intersections
int main()
{
    using ctype = double;
    using Box = Frackit::Box<ctype>;
    using Cylinder = Frackit::Cylinder<ctype>;
    using Circle = typename Cylinder::Circle;

    using Disk = Frackit::Disk<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;

    std::vector<ctype> scales({1.0e-3, 1.0, 1.0e3});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;
        const auto eps = Frackit::Precision<ctype>::confusion()*f;

        // make box with a cylindrical whole which we will use to create
        // parts of disks by computing the intersections
        const Box box(-10.0*f, -10.0*f, -10.0*f, 10.0*f, 10.0*f, 10.0*f);

        const Point center(0.0, 0.0, -f);
        const Direction axis(Vector(0.0, 0.0, 1.0));
        const Circle bottom(center, axis, f);
        const Cylinder cylinder(bottom, 2.0*f);

        const auto boxShape = Frackit::OCCUtilities::getShape(box);
        const auto cylinderShape = Frackit::OCCUtilities::getShape(cylinder);
        const auto domain = Frackit::OCCUtilities::cut(boxShape, cylinderShape, eps);

        // lambda to get the shapes of two disks confined to the domain
        auto getFaces = [&] (const auto& disk1, const auto& disk2) -> std::pair<TopoDS_Face, TopoDS_Face>
        {
            const auto diskShape1 = Frackit::OCCUtilities::getShape(disk1);
            const auto diskShape2 = Frackit::OCCUtilities::getShape(disk2);
            const auto diskCut1 = Frackit::OCCUtilities::intersect(diskShape1, domain, eps);
            const auto diskCut2 = Frackit::OCCUtilities::intersect(diskShape2, domain, eps);

            // get the faces of the result
            const auto diskCut1Faces = Frackit::OCCUtilities::getFaces(diskCut1);
            const auto diskCut2Faces = Frackit::OCCUtilities::getFaces(diskCut2);
            if (diskCut1Faces.size() != 1) throw std::runtime_error("Could not find single face of cut 1");
            if (diskCut2Faces.size() != 1) throw std::runtime_error("Could not find single face of cut 2");

            return {diskCut1Faces[0], diskCut2Faces[0]};
        };

        // create two disks and confine them to the domain to create new shapes
        const Point center1(0.0, -f, 0.0);
        const Direction majAxis1(Vector(0.0, 1.0, 0.0));
        const Direction minAxis1(Vector(1.0, 0.0, 0.0));
        const Disk disk1_1(center1, majAxis1, minAxis1, std::sqrt(2.0)*f, std::sqrt(2.0)*f);
        const Disk disk1_2(center1 + majAxis1*(2.0*f), majAxis1, minAxis1, std::sqrt(2.0)*f, std::sqrt(2.0)*f);

        auto [d1, d2] = getFaces(disk1_1, disk1_2);

        // intersection should give two points
        const Point expected1(f, 0.0, 0.0);
        const Point expected2(-f, 0.0, 0.0);
        auto result = Frackit::intersect(d1, d2, eps);
        if (result.size() != 2)
            throw std::runtime_error("Unexpected intersection size");
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result[0]);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result[1]);
        const auto& p1 = std::get<Point>(result[0]);
        const auto& p2 = std::get<Point>(result[1]);
        if (!p1.isEqual(expected1, eps) && !p1.isEqual(expected2, eps))
            throw std::runtime_error("Unexpected intersection point");
        if (!p2.isEqual(expected1, eps) && !p2.isEqual(expected2, eps))
            throw std::runtime_error("Unexpected intersection point");

        // create two faces such that the result gives two faces
        const Disk disk2_1(center1, majAxis1, minAxis1, 2.0*f, 2.0*f);
        const Disk disk2_2(center1 + majAxis1*(2.0*f), majAxis1, minAxis1, std::sqrt(2.0)*f, std::sqrt(2.0)*f);

        auto [d3, d4] = getFaces(disk2_1, disk2_2);
    BRepTools::Write(d3, "temp.brep");
    std::cout << "numwires1: " << Frackit::OCCUtilities::getWires(d3).size() << std::endl;
    std::cout << "numwires2: " << Frackit::OCCUtilities::getWires(d4).size() << std::endl;
        result = Frackit::intersect(d3, d4, eps);
        if (result.size() != 2)
            throw std::runtime_error("Unexpected intersection size");
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::face); }, result[0]);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::face); }, result[1]);

        // check face areas
        using std::abs;
        const auto deltaArea1 = abs(Frackit::computeArea(std::get<TopoDS_Face>(result[0])) - 0.26122571618784510851*f*f);
        std::cout << "deltaArea1 = " << deltaArea1 << ", eps = " << eps << std::endl;
        if (  deltaArea1 > eps )
            throw std::runtime_error("Unexpected intersection face area");

        const auto deltaArea2 = abs(Frackit::computeArea(std::get<TopoDS_Face>(result[1])) - 0.26122571618787740164*f*f);
        std::cout << "deltaArea2 = " << deltaArea2 << ", eps = " << eps << std::endl;
        if ( deltaArea2 > eps )
            throw std::runtime_error("Unexpected intersection face area");

        // create two faces such that the result gives two edges
        const Direction majAxis2(Vector(0.0, 1.0, 0.5));
        const Disk disk3_2(center1 - Vector(0.0, 0.0, 0.05*f), majAxis2, minAxis1, std::sqrt(2.0)*f, std::sqrt(2.0)*f);

        auto [d5, d6] = getFaces(disk2_1, disk3_2);
        result = Frackit::intersect(d5, d6, eps);
        if (result.size() != 2)
            throw std::runtime_error("Unexpected intersection size");
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::edge); }, result[0]);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::edge); }, result[1]);

        // check edge lengths
        const auto deltaLength1 = abs(Frackit::computeLength(std::get<TopoDS_Edge>(result[0])) - 0.97389732363386261049*f);
        std::cout << "deltaLength1 = " << deltaLength1 << ", eps = " << eps << std::endl;
        if (  deltaLength1 > eps )
            throw std::runtime_error("Unexpected intersection edge length");

        const auto deltaLength2 = abs(Frackit::computeLength(std::get<TopoDS_Edge>(result[1])) - 0.97389732363385483893*f);
        std::cout << "deltaLength2 = " << deltaLength2 << ", eps = " << eps << std::endl;
        if ( deltaLength2 > eps )
            throw std::runtime_error("Unexpected intersection edge length");

        // write result into .brep file for visualization
        TopoDS_Builder b;
        TopoDS_Compound c;
        b.MakeCompound(c);
        b.Add(c, d5);
        b.Add(c, d6);
        b.Add(c, std::get<TopoDS_Edge>(result[0]));
        b.Add(c, std::get<TopoDS_Edge>(result[1]));
        b.Add(c, domain);
        BRepTools::Write(c, "edge_intersections.brep");

        std::cout << "All tests passed"  << std::endl;
    }

    return 0;
}
