#include <vector>
#include <variant>
#include <type_traits>

#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopExp.hxx>

#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/polygon.hh>
#include <frackit/geometry/ellipsearc.hh>
#include <frackit/geometry/cylindersurface.hh>
#include <frackit/geometryutilities/name.hh>
#include <frackit/occ/breputilities.hh>

#include <frackit/magnitude/length.hh>
#include <frackit/intersection/intersect.hh>

enum IntersectionType { point, ellipseArc, ellipse, segment, empty };

template<class G>
void checkResultGeometry(const G& geometry, IntersectionType expected)
{
    throw std::runtime_error("Unexpected intersection geometry");
}

void checkResultGeometry(const Frackit::EmptyIntersection<3>& empty, IntersectionType expected)
{
    std::cout << "Found empty intersection" << std::endl;
    if (expected != IntersectionType::empty)
        throw std::runtime_error("Unexpected empty intersection");
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

template<class CT, int wd>
void checkResultGeometry(const Frackit::EllipseArc<CT, wd>& arc, IntersectionType expected)
{
    std::cout << "Found intersection ellipse arc with center " << arc.center() << ", "
              << "major axis: " << Frackit::Vector<CT, wd>(arc.majorAxis()) << ", "
              << "minor axis: " << Frackit::Vector<CT, wd>(arc.minorAxis()) << ", "
              << "major ax length: " << arc.majorAxisLength() << ", "
              << "minor ax length: " << arc.minorAxisLength() << ", "
              << "source: " << arc.source() << ", "
              << "target: " << arc.target() << ", "
              << "middle: " << arc.getPoint(0.5) << std::endl;
    if (expected != IntersectionType::ellipseArc)
        throw std::runtime_error("Got an unexpected arc intersection");
}

template<class CT, int wd>
void checkResultGeometry(const Frackit::Ellipse<CT, wd>& ellipse, IntersectionType expected)
{
    std::cout << "Found intersection ellipse with center " << ellipse.center() << ", "
              << "major axis: " << Frackit::Vector<CT, wd>(ellipse.majorAxis()) << ", "
              << "minor axis: " << Frackit::Vector<CT, wd>(ellipse.minorAxis()) << ", "
              << "major ax length: " << ellipse.majorAxisLength() << ", "
              << "minor ax length: " << ellipse.minorAxisLength() << std::endl;
    if (expected != IntersectionType::ellipse)
        throw std::runtime_error("Got an unexpected ellipse intersection");
}

template<class IS>
void checkResultGeometry(const std::vector<IS>& intersections, IntersectionType expected)
{
    for (const auto& is : intersections)
        checkResultGeometry(is, expected);
}

//! function to create the quadrilateral geometry from corner points
template<class QuadGeometry>
QuadGeometry makeQuad(const std::vector<Frackit::Point<typename QuadGeometry::ctype, 3>>& corners)
{
    using ctype = typename QuadGeometry::ctype;
    if constexpr (std::is_same_v<QuadGeometry, Frackit::Quadrilateral<ctype, 3>>)
        return QuadGeometry(corners[0], corners[1], corners[3], corners[2]);
    else if constexpr (std::is_same_v<QuadGeometry, Frackit::Polygon<ctype, 3>>)
        return QuadGeometry(corners);
    else
        throw std::runtime_error("Only Polygon or Quadrilateral are supported");
}

//! We test cylinder surface - quadrilateral intersections,
//! but also support quadrilaterals described by the polygon class
//! in order to test if the algorithm with polygon gives the same result
template<class QuadGeometry>
void doTest()
{
    using Quad = QuadGeometry;
    using ctype = typename Quad::ctype;
    using Point = typename Quad::Point;
    using CornerVector = std::vector<Point>;

    using Ellipse = Frackit::Ellipse<ctype, 3>;
    using Vector = Frackit::Vector<ctype, 3>;
    using EllipseArc = Frackit::EllipseArc<ctype, 3>;
    using CylinderSurface = Frackit::CylinderSurface<ctype>;
    using Segment = typename CylinderSurface::Segment;

    // OCC seems to fail below 1e-3 :(
    std::vector<ctype> scales({1.0e-3, 1.0, 1.0e3});
    for (auto f : scales)
    {
        std::cout << "\nChecking scale factor " << f << std::endl;

        Vector e1(1.0, 0.0, 0.0);
        Vector e2(0.0, 1.0, 0.0);
        Vector e3(0.0, 0.0, 1.0);
        CylinderSurface cylSurface(0.5*f, f);

        // epsilon for floating point comparison
        const ctype eps = Frackit::Precision<ctype>::confusion()*0.5*f;

        // quarilateral that touches the cylinder surface in one point
        auto quad1 = makeQuad<QuadGeometry>(CornerVector({Point(-0.5*f, 0.0*f, 0.0*f),
                                            Point(-1.0*f, 0.0*f, 0.0*f),
                                            Point(-1.0*f, 1.0*f, 0.0*f),
                                            Point(-0.5*f, 1.0*f, 0.0*f)}));
        auto result = intersect(cylSurface, quad1);
        auto resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad1));
        if (result.size() > 1)
            throw std::runtime_error("1: More than one intersection found");
        if (resultFromShape.size() > 1)
            throw std::runtime_error("1: More than one intersection found (shapes)");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, variant);
        for (const auto& variant : resultFromShape)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, variant);
        const auto p1 = std::get<Point>(result[0]);
        if (!p1.isEqual(Point({-0.5*f, 0.0, 0.0}), eps))
            throw std::runtime_error("Unexpected point 1");
        if (!p1.isEqual(std::get<Point>(resultFromShape[0]), eps))
            throw std::runtime_error("Point obtained from shapes not equal to that from internal geoms");
        std::cout << "Test passed" << std::endl;

        // quad that describes a horizontal cross-section (circle) of the surface
        auto quad2 = makeQuad<QuadGeometry>(CornerVector({Point(-1.0*f, -1.0*f, 0.0*f),
                                            Point( 1.0*f, -1.0*f, 0.0*f),
                                            Point( 1.0*f, 1.0*f, 0.0*f),
                                            Point(-1.0*f, 1.0*f, 0.0*f)}));
        result = intersect(cylSurface, quad2);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad2));
        if (result.size() > 1)
            throw std::runtime_error("2: More than one intersection found");
        if (resultFromShape.size() > 1)
            throw std::runtime_error("2: More than one intersection found (shapes)");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipse); }, variant);
        if (!std::holds_alternative<TopoDS_Edge>(resultFromShape[0]))
            throw std::runtime_error("Unexpected result type for is from shapes");
        const auto ellipse1 = std::get<Ellipse>(result[0]);
        using std::abs;
        if (abs(Vector(ellipse1.normal())*e1) > eps)
            throw std::runtime_error("Unexpected normal direction of ellipse 1");
        if (abs(Vector(ellipse1.normal())*e2) > eps)
            throw std::runtime_error("Unexpected normal direction of ellipse 1");
        if (abs(1.0 - abs(Vector(ellipse1.normal())*e3)) > eps)
            throw std::runtime_error("Unexpected normal direction of ellipse 1");

        TopoDS_Vertex v1, v2;
        TopExp::Vertices(std::get<TopoDS_Edge>(resultFromShape[0]), v1, v2);
        auto pShape1 = Frackit::OCCUtilities::point(v1);
        auto pShape2 = Frackit::OCCUtilities::point(v2);
        if (!pShape1.isEqual(pShape2, eps))
            throw std::runtime_error("Shape intersection did not yield ellipse");
        if (!ellipse1.contains(pShape1, eps) || !ellipse1.contains(pShape2, eps))
            throw std::runtime_error("Shape intersection points not on ellipse!");
        std::cout << "Test passed" << std::endl;

        // quad that describes a vertical cross-section of the surface (two segments)
        auto quad3 = makeQuad<QuadGeometry>(CornerVector({Point(-1.0*f, 0.0*f, -1.0*f),
                                            Point(-1.0*f, 0.0*f,  1.0*f),
                                            Point( 1.0*f, 0.0*f,  1.0*f),
                                            Point( 1.0*f, 0.0*f, -1.0*f)}));
        result = intersect(cylSurface, quad3);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad3));
        if (result.size() != 2)
            throw std::runtime_error("3: Did not find two intersection geometries");
        if (resultFromShape.size() != 2)
            throw std::runtime_error("3: Did not find two intersection geometries (shapes)");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, variant);
        for (const auto& variant : resultFromShape)
            if (!std::holds_alternative<TopoDS_Edge>(variant))
                throw std::runtime_error("Did not obtain 2 edges from intersection with shapes");
        if ( abs(Frackit::computeLength(std::get<Segment>(result[0])) - cylSurface.height()) > eps )
            throw std::runtime_error("Unexpected segment height (test 3)");
        if ( abs(Frackit::computeLength(std::get<Segment>(result[1])) - cylSurface.height()) > eps )
            throw std::runtime_error("Unexpected segment height (test 3)");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Segment>(is).getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected segment (test 3)");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Segment>(is).getPoint(0.5).isEqual(Point(-0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected segment (test 3)");

        TopExp::Vertices(std::get<TopoDS_Edge>(resultFromShape[0]), v1, v2);
        pShape1 = Frackit::OCCUtilities::point(v1);
        pShape2 = Frackit::OCCUtilities::point(v2);
        if (abs(Segment(pShape1, pShape2).length() - cylSurface.height()) > eps)
            throw std::runtime_error("Edge length from shape intersections has incorrect length");
        TopExp::Vertices(std::get<TopoDS_Edge>(resultFromShape[1]), v1, v2);
        pShape1 = Frackit::OCCUtilities::point(v1);
        pShape2 = Frackit::OCCUtilities::point(v2);
        if (abs(Segment(pShape1, pShape2).length() - cylSurface.height()) > eps)
            throw std::runtime_error("Edge length from shape intersections has incorrect length");

        std::cout << "Test passed" << std::endl;

        // quad that touches the cylinder surface in one point (opposite side)
        auto quad4 = makeQuad<QuadGeometry>(CornerVector({Point(0.5*f, 0.0*f, 0.0*f),
                                            Point(1.0*f, 0.0*f, 0.0*f),
                                            Point(1.0*f, 1.0*f, 0.0*f),
                                            Point(0.5*f, 1.0*f, 0.0*f)}));
        result = intersect(cylSurface, quad4);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad4));
        if (result.size() > 1)
            throw std::runtime_error("4: More than one intersection found");
        if (resultFromShape.size() > 1)
            throw std::runtime_error("4: More than one intersection found (shapes)");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, variant);
        for (const auto& variant : resultFromShape)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, variant);
        if (!std::get<Point>(result[0]).isEqual(Point({0.5*f, 0.0, 0.0}), eps))
            throw std::runtime_error("Unexpected point (test 4)");
        if (!std::get<Point>(resultFromShape[0]).isEqual(Point({0.5*f, 0.0, 0.0}), eps))
            throw std::runtime_error("Unexpected point (test 4 - shapes)");
        std::cout << "Test passed" << std::endl;

        // no intersection
        auto quad5 = makeQuad<QuadGeometry>(CornerVector({Point(0.5*f+1e-3*f, 0.0*f, 0.0*f),
                                            Point(1.0*f, 0.0*f, 0.0*f),
                                            Point(1.0*f, 1.0*f, 0.0*f),
                                            Point(0.5*f, 1.0*f, 0.0*f)}));
        result = intersect(cylSurface, quad5);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad5));
        if (result.size() > 1)
            throw std::runtime_error("5: More than one intersection found");
        if (resultFromShape.size() > 1)
            throw std::runtime_error("5: More than one intersection found");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::empty); }, variant);
        for (const auto& variant : resultFromShape)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::empty); }, variant);
        std::cout << "Test passed" << std::endl;

        // quad that intersects on one side only
        auto quad6 = makeQuad<QuadGeometry>(CornerVector({Point(0.0*f,  0.0*f, 0.5*f),
                                            Point(1.0*f, -1.0*f, 0.5*f),
                                            Point(2.0*f,  0.0*f, 0.5*f),
                                            Point(1.0*f,  1.0*f, 0.5*f)}));
        result = intersect(cylSurface, quad6);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad6));
        if (result.size() > 1)
            throw std::runtime_error("6: More than one intersection found");
        if (resultFromShape.size() != 2) // gives two edges as it intersects with closing segment in shape representation
            throw std::runtime_error("6: More than one intersection found (shapes)");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, variant);
        if (!std::holds_alternative<TopoDS_Edge>(resultFromShape[0]))
            throw std::runtime_error("Expected edge intersection in test 6");
        if (!std::holds_alternative<TopoDS_Edge>(resultFromShape[1]))
            throw std::runtime_error("Expected edge intersection in test 6");
        const auto arc6 = std::get<EllipseArc>(result[0]);
        if (!arc6.getPoint(0.5).isEqual(Point({0.5*f, 0.0, 0.5*f}), eps))
            throw std::runtime_error("Unexpected ellipse arc (test 6)");

        TopExp::Vertices(std::get<TopoDS_Edge>(resultFromShape[0]), v1, v2);
        pShape1 = Frackit::OCCUtilities::point(v1);
        pShape2 = Frackit::OCCUtilities::point(v2);
        if (!pShape1.isEqual(arc6.source(), eps) && !pShape1.isEqual(arc6.target(), eps)
            && !pShape2.isEqual(arc6.source(), eps) && !pShape2.isEqual(arc6.target(), eps))
            throw std::runtime_error("Point is not on ellipse arc of test 6");

        TopExp::Vertices(std::get<TopoDS_Edge>(resultFromShape[1]), v1, v2);
        pShape1 = Frackit::OCCUtilities::point(v1);
        pShape2 = Frackit::OCCUtilities::point(v2);
        if (!pShape1.isEqual(arc6.source(), eps) && !pShape1.isEqual(arc6.target(), eps)
            && !pShape2.isEqual(arc6.source(), eps) && !pShape2.isEqual(arc6.target(), eps))
            throw std::runtime_error("Point is not on ellipse arc of test 6");
        std::cout << "Test passed" << std::endl;

        // quad that intersects on two sides
        auto quad7 = makeQuad<QuadGeometry>(CornerVector({Point(-2.0*f, -0.25*f, 0.5*f),
                                            Point( 2.0*f, -0.25*f, 0.5*f),
                                            Point( 2.0*f,  0.25*f, 0.5*f),
                                            Point(-2.0*f,  0.25*f, 0.5*f)}));
        result = intersect(cylSurface, quad7);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad7));
        if (result.size() != 2)
            throw std::runtime_error("7: Did not find two intersections");
        if (resultFromShape.size() != 3) // gives three edges as it intersects with closing segment in shape representation
            throw std::runtime_error("7: Did not find two intersections");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, variant);
        for (const auto& variant : resultFromShape)
            if (!std::holds_alternative<TopoDS_Edge>(variant))
                throw std::runtime_error("Expected edge intersections in test 7");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected ellipse arc (test 7)");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(-0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected ellipse arc (test 7)");

        const auto arc7_1 = std::get<EllipseArc>(result[0]);
        const auto arc7_2 = std::get<EllipseArc>(result[1]);
        for (const auto& shapeResult : resultFromShape)
        {
            TopExp::Vertices(std::get<TopoDS_Edge>(shapeResult), v1, v2);
            pShape1 = Frackit::OCCUtilities::point(v1);
            pShape2 = Frackit::OCCUtilities::point(v2);
            if (!pShape1.isEqual(arc7_1.source(), eps) && !pShape1.isEqual(arc7_1.target(), eps)
                && !pShape1.isEqual(arc7_2.source(), eps) && !pShape1.isEqual(arc7_2.target(), eps)
                && !pShape2.isEqual(arc7_1.source(), eps) && !pShape2.isEqual(arc7_1.target(), eps)
                && !pShape2.isEqual(arc7_2.source(), eps) && !pShape2.isEqual(arc7_2.target(), eps))
                throw std::runtime_error("Point is not on ellipse arc of test 7");
        }

        std::cout << "Test passed" << std::endl;

        // same as the previous test but with an inclined quad
        auto quad8 = makeQuad<QuadGeometry>(CornerVector({Point(-2.0*f, -0.25*f, 0.5*f + 0.3*2.0*f),
                                            Point( 2.0*f, -0.25*f, 0.5*f - 0.3*2.0*f),
                                            Point( 2.0*f,  0.25*f, 0.5*f - 0.3*2.0*f),
                                            Point(-2.0*f,  0.25*f, 0.5*f + 0.3*2.0*f)}));
        result = intersect(cylSurface, quad8);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad8));
        if (result.size() != 2)
            throw std::runtime_error("8: Did not find two intersections");
        if (resultFromShape.size() != 3) // gives three edges as it intersects with closing segment in shape representation
            throw std::runtime_error("8: Did not find two intersections");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, variant);
        for (const auto& variant : resultFromShape)
            if (!std::holds_alternative<TopoDS_Edge>(variant))
                throw std::runtime_error("Expected edge shapes in test 8");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point( 0.5*f, 0.0, 0.5*f - 0.3*0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected ellipse arc (test 8)");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(-0.5*f, 0.0, 0.5*f + 0.3*0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected ellipse arc (test 8)");

        const auto arc8_1 = std::get<EllipseArc>(result[0]);
        const auto arc8_2 = std::get<EllipseArc>(result[1]);
        for (const auto& shapeResult : resultFromShape)
        {
            TopExp::Vertices(std::get<TopoDS_Edge>(shapeResult), v1, v2);
            pShape1 = Frackit::OCCUtilities::point(v1);
            pShape2 = Frackit::OCCUtilities::point(v2);
            if (!pShape1.isEqual(arc8_1.source(), eps) && !pShape1.isEqual(arc8_1.target(), eps)
                && !pShape1.isEqual(arc8_2.source(), eps) && !pShape1.isEqual(arc8_2.target(), eps)
                && !pShape2.isEqual(arc8_1.source(), eps) && !pShape2.isEqual(arc8_1.target(), eps)
                && !pShape2.isEqual(arc8_2.source(), eps) && !pShape2.isEqual(arc8_2.target(), eps))
                throw std::runtime_error("Point is not on ellipse arc of test 8");
        }

        std::cout << "Test passed" << std::endl;

        // quad that is much bigger, intersection an ellipse
        auto quad9 = makeQuad<QuadGeometry>(CornerVector({Point(-2.0*f, -2.0*f, 0.5*f + 0.3*2.0*f),
                                            Point( 2.0*f, -2.0*f, 0.5*f - 0.3*2.0*f),
                                            Point( 2.0*f,  2.0*f, 0.5*f - 0.3*2.0*f),
                                            Point(-2.0*f,  2.0*f, 0.5*f + 0.3*2.0*f)}));
        result = intersect(cylSurface, quad9);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad9));
        if (result.size() > 1)
            throw std::runtime_error("9: Found more than one intersection");
        if (resultFromShape.size() > 1)
            throw std::runtime_error("9: Found more than one intersection (shapes)");
        std::cout << "Checking test 9" << std::endl;
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipse); }, variant);
        if (!std::get<Ellipse>(result[0]).getPoint(0.0).isEqual(Point({-0.5*f, 0.0, 0.5*f + 0.3*0.5*f}), eps))
            throw std::runtime_error("Unexpected ellipse (test 9)");
        if (!std::get<Ellipse>(result[0]).getPoint(0.5).isEqual(Point({0.5*f, 0.0, 0.5*f - 0.3*0.5*f}), eps))
            throw std::runtime_error("Unexpected ellipse (test 9)");
        if (!std::holds_alternative<TopoDS_Edge>(resultFromShape[0]))
            throw std::runtime_error("Expected edge shape as intersection");

        const auto ellipse_8 = std::get<Ellipse>(result[0]);
        for (const auto& shapeResult : resultFromShape)
        {
            TopExp::Vertices(std::get<TopoDS_Edge>(shapeResult), v1, v2);
            pShape1 = Frackit::OCCUtilities::point(v1);
            pShape2 = Frackit::OCCUtilities::point(v2);
            if (!ellipse_8.contains(pShape1, eps) || !ellipse_8.contains(pShape2, eps))
                throw std::runtime_error("Points from shape intersections not on ellipse!");
        }

        std::cout << "Test passed" << std::endl;

        // quad that touches in two points
        auto quad10 = makeQuad<QuadGeometry>(CornerVector({Point(-0.5*f,  0.0*f, 0.5*f),
                                             Point( 0.0*f, -0.25*f, 0.5*f),
                                             Point( 0.5*f,  0.0*f, 0.5*f),
                                             Point( 0.0*f,  0.25*f, 0.5*f)}));
        result = intersect(cylSurface, quad10);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad10));
        if (result.size() != 2)
            throw std::runtime_error("10: Did not find two intersections");
        if (resultFromShape.size() != 2)
            throw std::runtime_error("10: Did not find two intersections");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, variant);
        for (const auto& variant : resultFromShape)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, variant);
        if (!std::get<Point>(result[0]).isEqual(Point({-0.5*f, 0.0, 0.5*f}), eps))
            throw std::runtime_error("Unexpected point 10");
        if (!std::get<Point>(result[1]).isEqual(Point({0.5*f, 0.0, 0.5*f}), eps))
            throw std::runtime_error("Unexpected point 10");

        pShape1 = std::get<Point>(resultFromShape[0]);
        pShape2 = std::get<Point>(resultFromShape[1]);
        if (!pShape1.isEqual(std::get<Point>(result[0]), eps)
            && !pShape1.isEqual(std::get<Point>(result[1]), eps))
            throw std::runtime_error("Wrong point from shape intersection");
        if (!pShape2.isEqual(std::get<Point>(result[0]), eps)
            && !pShape2.isEqual(std::get<Point>(result[1]), eps))
            throw std::runtime_error("Wrong point from shape intersection");

        std::cout << "Test passed" << std::endl;

        // quad that intersects in four points
        auto quad11 = makeQuad<QuadGeometry>(CornerVector({Point(-0.5*f,  0.0*f, 0.5*f),
                                             Point( 0.0*f, -0.5*f, 0.5*f),
                                             Point( 0.5*f,  0.0*f, 0.5*f),
                                             Point( 0.0*f,  0.5*f, 0.5*f)}));
        result = intersect(cylSurface, quad11);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad11));
        if (result.size() != 4)
            throw std::runtime_error("11: Did not find four intersections");
        if (resultFromShape.size() != 4)
            throw std::runtime_error("11: Did not find four intersections");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, variant);
        for (const auto& variant : resultFromShape)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, variant);
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Point>(is).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected point (test 11)");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Point>(is).isEqual(Point(0.0*f, 0.5*f, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected point (test 11)");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Point>(is).isEqual(Point(-0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected point (test 11)");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Point>(is).isEqual(Point(0.0*f, -0.5*f, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected point (test 11)");

        if (std::count_if(resultFromShape.begin(),
                          resultFromShape.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Point>(is).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected point (test 11 - shapes)");
        if (std::count_if(resultFromShape.begin(),
                          resultFromShape.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Point>(is).isEqual(Point(0.0*f, 0.5*f, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected point (test 11 - shapes)");
        if (std::count_if(resultFromShape.begin(),
                          resultFromShape.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Point>(is).isEqual(Point(-0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected point (test 11 - shapes)");
        if (std::count_if(resultFromShape.begin(),
                          resultFromShape.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Point>(is).isEqual(Point(0.0*f, -0.5*f, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected point (test 11 - shapes)");

        std::cout << "Test passed" << std::endl;

        // quad that intersects in two short segments
        auto quad12 = makeQuad<QuadGeometry>(CornerVector({Point(-0.6*f, 0.0*f, 0.5*f),
                                             Point( 0.0*f, 0.0*f, 0.4*f),
                                             Point( 0.6*f, 0.0*f, 0.5*f),
                                             Point( 0.0*f, 0.0*f, 0.6*f)}));
        result = intersect(cylSurface, quad12);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad12));
        if (result.size() != 2)
            throw std::runtime_error("12: Did not find two intersections");
        if (resultFromShape.size() != 2)
            throw std::runtime_error("12: Did not find two intersections");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, variant);
        for (const auto& variant : resultFromShape)
            if (!std::holds_alternative<TopoDS_Edge>(variant))
                throw std::runtime_error("Expected edge intersections in test 12");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Segment>(is).getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected segment center (test 12)");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Segment>(is).getPoint(0.5).isEqual(Point(-0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected segment center (test 12)");

        for (const auto& shapeResult : resultFromShape)
        {
            TopExp::Vertices(std::get<TopoDS_Edge>(shapeResult), v1, v2);
            pShape1 = Frackit::OCCUtilities::point(v1);
            pShape2 = Frackit::OCCUtilities::point(v2);
            auto tmp = Segment(pShape1, pShape2).getPoint(0.5);
            if (!tmp.isEqual(Point(0.5*f, 0.0, 0.5*f), eps)
                && !tmp.isEqual(Point(-0.5*f, 0.0, 0.5*f), eps))
                throw std::runtime_error("Shape intersection segment incorrect (test 12)");
        }

        std::cout << "Test passed" << std::endl;

        // quad that intersects in two short arcs
        auto quad13 = makeQuad<QuadGeometry>(CornerVector({Point(-0.6*f,  0.0*f, 0.5*f),
                                             Point( 0.0*f, -0.1*f, 0.5*f),
                                             Point( 0.6*f,  0.0*f, 0.5*f),
                                             Point( 0.0*f,  0.1*f, 0.5*f)}));
        result = intersect(cylSurface, quad13);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad13));
        if (result.size() != 2)
            throw std::runtime_error("13: Did not find two intersections");
        if (resultFromShape.size() != 3) // gives three edges as it intersects with closing segment in shape representation
            throw std::runtime_error("13: Did not find two intersections (shape)");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, variant);
        for (const auto& variant : resultFromShape)
            if (!std::holds_alternative<TopoDS_Edge>(variant))
                throw std::runtime_error("Expected edge intersections (test 13)");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected ellipse arc (test 13)");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(-0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected ellipse arc (test 13)");

        for (const auto& shapeResult : resultFromShape)
        {
            TopExp::Vertices(std::get<TopoDS_Edge>(shapeResult), v1, v2);
            pShape1 = Frackit::OCCUtilities::point(v1);
            pShape2 = Frackit::OCCUtilities::point(v2);
            std::vector<Point> arcTips; arcTips.reserve(4);
            arcTips.emplace_back(std::get<EllipseArc>(result[0]).source());
            arcTips.emplace_back(std::get<EllipseArc>(result[1]).target());
            arcTips.emplace_back(std::get<EllipseArc>(result[1]).source());
            arcTips.emplace_back(std::get<EllipseArc>(result[0]).target());
            if (std::none_of(arcTips.begin(), arcTips.end(),
                             [&] (const auto& tip) { return tip.isEqual(pShape1, eps); })
                && std::none_of(arcTips.begin(), arcTips.end(),
                                [&] (const auto& tip) { return tip.isEqual(pShape2, eps); }))
                throw std::runtime_error("Shape intersection point incorrect (test 13)");
        }

        std::cout << "Test passed" << std::endl;

        // quad that intersects in two long arcs
        auto quad14 = makeQuad<QuadGeometry>(CornerVector({Point(-0.6*f,  0.0*f, 0.5*f),
                                             Point( 0.0*f, -0.4*f, 0.5*f),
                                             Point( 0.6*f,  0.0*f, 0.5*f),
                                             Point( 0.0*f,  0.4*f, 0.5*f)}));
        result = intersect(cylSurface, quad14);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad14));
        if (result.size() != 2)
            throw std::runtime_error("14: Did not find two intersections");
        if (resultFromShape.size() != 3) // gives three edges as it intersects with closing segment in shape representation
            throw std::runtime_error("14: Did not find two intersections (shapes)");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, variant);
        for (const auto& variant : resultFromShape)
            if (!std::holds_alternative<TopoDS_Edge>(variant))
                throw std::runtime_error("Expected edge shape as intersection");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected ellipse arc (test 14)");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected ellipse arc (test 14)");

        for (const auto& shapeResult : resultFromShape)
        {
            TopExp::Vertices(std::get<TopoDS_Edge>(shapeResult), v1, v2);
            pShape1 = Frackit::OCCUtilities::point(v1);
            pShape2 = Frackit::OCCUtilities::point(v2);
            std::vector<Point> arcTips; arcTips.reserve(4);
            arcTips.emplace_back(std::get<EllipseArc>(result[0]).source());
            arcTips.emplace_back(std::get<EllipseArc>(result[1]).target());
            arcTips.emplace_back(std::get<EllipseArc>(result[1]).source());
            arcTips.emplace_back(std::get<EllipseArc>(result[0]).target());
            if (std::none_of(arcTips.begin(), arcTips.end(),
                             [&] (const auto& tip) { return tip.isEqual(pShape1, eps); })
                && std::none_of(arcTips.begin(), arcTips.end(),
                                [&] (const auto& tip) { return tip.isEqual(pShape2, eps); }))
                throw std::runtime_error("Shape intersection point incorrect (test 14)");
        }

        std::cout << "Test passed" << std::endl;

        // quad that intersects in one segment and one point
        auto quad15 = makeQuad<QuadGeometry>(CornerVector({Point(-0.6*f, 0.0*f, 0.5*f),
                                             Point( 0.0*f, 0.0*f, 0.4*f),
                                             Point( 0.5*f, 0.0*f, 0.5*f),
                                             Point( 0.0*f, 0.0*f, 0.6*f)}));
        result = intersect(cylSurface, quad15);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad15));
        if (result.size() != 2)
            throw std::runtime_error("15: Did not find two intersections");
        if (resultFromShape.size() != 2)
            throw std::runtime_error("15: Did not find two intersections");
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result[0]);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result[1]);
        if (!std::holds_alternative<Point>(resultFromShape[0])
            && !std::holds_alternative<TopoDS_Edge>(resultFromShape[0]))
            throw std::runtime_error("Unexpected intersection shape (test 15)");
        if (!std::holds_alternative<Point>(resultFromShape[1])
            && !std::holds_alternative<TopoDS_Edge>(resultFromShape[1]))
            throw std::runtime_error("Unexpected intersection shape (test 15)");

        std::cout << "Test passed" << std::endl;

        // quad that intersects in one arc and one point
        auto quad16 = makeQuad<QuadGeometry>(CornerVector({Point(-0.6*f,  0.0*f, 0.5*f),
                                             Point( 0.0*f, -0.1*f, 0.5*f),
                                             Point( 0.5*f,  0.0*f, 0.5*f),
                                             Point( 0.0*f,  0.1*f, 0.5*f)}));
        result = intersect(cylSurface, quad16);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad16));
        if (result.size() != 2)
            throw std::runtime_error("16: Did not find two intersections");
        if (resultFromShape.size() != 2)
            throw std::runtime_error("16: Did not find two intersections (from shapes)");
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result[0]);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, result[1]);
        if (!std::holds_alternative<Point>(resultFromShape[0])
            && !std::holds_alternative<TopoDS_Edge>(resultFromShape[0]))
            throw std::runtime_error("Unexpected intersection shape (test 16)");
        if (!std::holds_alternative<Point>(resultFromShape[1])
            && !std::holds_alternative<TopoDS_Edge>(resultFromShape[1]))
            throw std::runtime_error("Unexpected intersection shape (test 16)");
        std::cout << "Test passed" << std::endl;

        // disk that intersects in one parabola
        auto quad17 = makeQuad<QuadGeometry>(CornerVector({Point( 0.0*f, -1.0*f, 0.0*f),
                                             Point( 2.0*f, -0.5*f, 0.5*f),
                                             Point( 0.0*f,  1.0*f, 2.0*f),
                                             Point(-2.0*f, -0.5*f, 0.5*f)}));
        result = intersect(cylSurface, quad17);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad17));
        if (result.size() != 1)
            throw std::runtime_error("17: Did not find a single intersections");
        if (resultFromShape.size() != 1)
            throw std::runtime_error("17: Did not find a single intersections");
        if (abs(std::get<EllipseArc>(result[0]).source().z() - f) > eps)
            throw std::runtime_error("Unexpected ellipse arc (source point)");
        if (abs(std::get<EllipseArc>(result[0]).target().z() - f) > eps)
            throw std::runtime_error("Unexpected ellipse arc (target point)");
        if (abs(std::get<EllipseArc>(result[0]).getPoint(0.5).x()) > eps)
            throw std::runtime_error("Unexpected ellipse arc");
        if (abs(std::get<EllipseArc>(result[0]).getPoint(0.5).y() + 0.5*f) > eps)
            throw std::runtime_error("Unexpected ellipse arc");
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, result[0]);

        if (!std::holds_alternative<TopoDS_Edge>(resultFromShape[0]))
            throw std::runtime_error("Unexpected intersection shape (test 17)");

        for (const auto& shapeResult : resultFromShape)
        {
            TopExp::Vertices(std::get<TopoDS_Edge>(shapeResult), v1, v2);
            pShape1 = Frackit::OCCUtilities::point(v1);
            pShape2 = Frackit::OCCUtilities::point(v2);
            std::vector<Point> arcTips; arcTips.reserve(4);
            arcTips.emplace_back(std::get<EllipseArc>(result[0]).source());
            arcTips.emplace_back(std::get<EllipseArc>(result[0]).target());
            if (std::none_of(arcTips.begin(), arcTips.end(),
                             [&] (const auto& tip) { return tip.isEqual(pShape1, eps); })
                && std::none_of(arcTips.begin(), arcTips.end(),
                                [&] (const auto& tip) { return tip.isEqual(pShape2, eps); }))
                throw std::runtime_error("Shape intersection point incorrect (test 17)");
        }

        std::cout << "Test passed" << std::endl;

        // quad that touches in a segment
        auto quad18 = makeQuad<QuadGeometry>(CornerVector({Point(0.0*f, 0.5*f, 0.25*f),
                                             Point(0.0*f, 0.5*f, 0.75*f),
                                             Point(0.0*f, 1.0*f, 0.75*f),
                                             Point(0.0*f, 1.0*f, 0.25*f)}));
        result = intersect(cylSurface, quad18);
        resultFromShape = intersect(cylSurface, Frackit::OCCUtilities::getShape(quad18));
        if (result.size() != 1)
            throw std::runtime_error("18: Did not find a single intersections");
        if (resultFromShape.size() != 1)
            throw std::runtime_error("18: Did not find a single intersections (shapes)");
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result[0]);
        if (abs(std::get<Segment>(result[0]).source().x() - 0.0) > eps)
            throw std::runtime_error("Unexpected segment (source point)");
        if (abs(std::get<Segment>(result[0]).target().x() - 0.0) > eps)
            throw std::runtime_error("Unexpected segment (target point)");
        if (abs(std::get<Segment>(result[0]).getPoint(0.5).z() - 0.5*f) > eps)
            throw std::runtime_error("Unexpected segment (mid point)");

        if (!std::holds_alternative<TopoDS_Edge>(resultFromShape[0]))
            throw std::runtime_error("Unexpected intersection shape (test 18)");

        for (const auto& shapeResult : resultFromShape)
        {
            TopExp::Vertices(std::get<TopoDS_Edge>(shapeResult), v1, v2);
            pShape1 = Frackit::OCCUtilities::point(v1);
            pShape2 = Frackit::OCCUtilities::point(v2);
            std::vector<Point> segmentTips; segmentTips.reserve(4);
            segmentTips.emplace_back(std::get<Segment>(result[0]).source());
            segmentTips.emplace_back(std::get<Segment>(result[0]).target());
            if (std::none_of(segmentTips.begin(), segmentTips.end(),
                             [&] (const auto& tip) { return tip.isEqual(pShape1, eps); })
                && std::none_of(segmentTips.begin(), segmentTips.end(),
                                [&] (const auto& tip) { return tip.isEqual(pShape2, eps); }))
                throw std::runtime_error("Shape intersection point incorrect (test 18)");
        }

        std::cout << "Test passed" << std::endl;
        std::cout << "All tests passed"  << std::endl;
    }
}

//! test both with quadrilaterals and polygons
int main()
{
    using ctype = double;
    using Quad = Frackit::Quadrilateral<ctype, 3>;
    using Polygon = Frackit::Polygon<ctype, 3>;

    doTest<Quad>();
    doTest<Polygon>();

    return 0;
}
