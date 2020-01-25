#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/ellipsearc.hh>
#include <frackit/geometry/cylindersurface.hh>

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

//! test cylinder surface - quadrilateral intersections
int main()
{
    using ctype = double;
    using CylinderSurface = Frackit::CylinderSurface<ctype>;
    using Segment = typename CylinderSurface::Segment;
    using Quad = Frackit::Quadrilateral<ctype, 3>;
    using Point = typename Quad::Point;
    using Ellipse = Frackit::Ellipse<ctype, 3>;
    using Vector = Frackit::Vector<ctype, 3>;
    using EllipseArc = Frackit::EllipseArc<ctype, 3>;

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
        Quad quad1(Point(-0.5*f, 0.0*f, 0.0*f),
                   Point(-1.0*f, 0.0*f, 0.0*f),
                   Point(-0.5*f, 1.0*f, 0.0*f),
                   Point(-1.0*f, 1.0*f, 0.0*f));
        auto result = intersect(cylSurface, quad1);
        if (result.size() > 1)
            throw std::runtime_error("1: More than one intersection found");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, variant);
        const auto p1 = std::get<Point>(result[0]);
        if (!p1.isEqual(Point({-0.5*f, 0.0, 0.0}), eps))
            throw std::runtime_error("Unexpected point 1");
        std::cout << "Test passed" << std::endl;

        // quad that describes a horizontal cross-section (circle) of the surface
        Quad quad2(Point(-1.0*f, -1.0*f, 0.0*f),
                   Point( 1.0*f, -1.0*f, 0.0*f),
                   Point(-1.0*f, 1.0*f, 0.0*f),
                   Point( 1.0*f, 1.0*f, 0.0*f));
        result = intersect(cylSurface, quad2);
        if (result.size() > 1)
            throw std::runtime_error("2: More than one intersection found");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipse); }, variant);
        const auto ellipse1 = std::get<Ellipse>(result[0]);
        using std::abs;
        if (abs(Vector(ellipse1.normal())*e1) > eps)
            throw std::runtime_error("Unexpected normal direction of ellipse 1");
        if (abs(Vector(ellipse1.normal())*e2) > eps)
            throw std::runtime_error("Unexpected normal direction of ellipse 1");
        if (abs(1.0 - abs(Vector(ellipse1.normal())*e3)) > eps)
            throw std::runtime_error("Unexpected normal direction of ellipse 1");
        std::cout << "Test passed" << std::endl;

        // quad that describes a vertical cross-section of the surface (two segments)
        Quad quad3(Point(-1.0*f, 0.0*f, -1.0*f),
                   Point(-1.0*f, 0.0*f,  1.0*f),
                   Point( 1.0*f, 0.0*f, -1.0*f),
                   Point( 1.0*f, 0.0*f,  1.0*f));
        result = intersect(cylSurface, quad3);
        if (result.size() != 2)
            throw std::runtime_error("3: Did not find two intersection geometries");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, variant);
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
        std::cout << "Test passed" << std::endl;

        // quad that touches the cylinder surface in one point (opposite side)
        Quad quad4(Point(0.5*f, 0.0*f, 0.0*f),
                   Point(1.0*f, 0.0*f, 0.0*f),
                   Point(0.5*f, 1.0*f, 0.0*f),
                   Point(1.0*f, 1.0*f, 0.0*f));
        result = intersect(cylSurface, quad4);
        if (result.size() > 1)
            throw std::runtime_error("4: More than one intersection found");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, variant);
        if (!std::get<Point>(result[0]).isEqual(Point({0.5*f, 0.0, 0.0}), eps))
            throw std::runtime_error("Unexpected point (test 4)");
        std::cout << "Test passed" << std::endl;

        // no intersection
        Quad quad5(Point(0.5*f+1e-3*f, 0.0*f, 0.0*f),
                   Point(1.0*f, 0.0*f, 0.0*f),
                   Point(0.5*f, 1.0*f, 0.0*f),
                   Point(1.0*f, 1.0*f, 0.0*f));
        result = intersect(cylSurface, quad5);
        if (result.size() > 1)
            throw std::runtime_error("5: More than one intersection found");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::empty); }, variant);
        std::cout << "Test passed" << std::endl;

        // quad that intersects on one side only
        Quad quad6(Point(0.0*f,  0.0*f, 0.5*f),
                   Point(1.0*f, -1.0*f, 0.5*f),
                   Point(1.0*f,  1.0*f, 0.5*f),
                   Point(2.0*f,  0.0*f, 0.5*f));
        result = intersect(cylSurface, quad6);
        if (result.size() > 1)
            throw std::runtime_error("6: More than one intersection found");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, variant);
        if (!std::get<EllipseArc>(result[0]).getPoint(0.5).isEqual(Point({0.5*f, 0.0, 0.5*f}), eps))
            throw std::runtime_error("Unexpected ellipse arc (test 6)");
        std::cout << "Test passed" << std::endl;

        // quad that intersects on two sides
        Quad quad7(Point(-2.0*f, -0.25*f, 0.5*f),
                   Point( 2.0*f, -0.25*f, 0.5*f),
                   Point(-2.0*f,  0.25*f, 0.5*f),
                   Point( 2.0*f,  0.25*f, 0.5*f));
        result = intersect(cylSurface, quad7);
        if (result.size() != 2)
            throw std::runtime_error("7: Did not find two intersections");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, variant);
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
        std::cout << "Test passed" << std::endl;

        // same as the previous test but with an inclined quad
        Quad quad8(Point(-2.0*f, -0.25*f, 0.5*f + 0.3*2.0*f),
                   Point( 2.0*f, -0.25*f, 0.5*f - 0.3*2.0*f),
                   Point(-2.0*f,  0.25*f, 0.5*f + 0.3*2.0*f),
                   Point( 2.0*f,  0.25*f, 0.5*f - 0.3*2.0*f));
        result = intersect(cylSurface, quad8);
        if (result.size() != 2)
            throw std::runtime_error("8: Did not find two intersections");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, variant);
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
        std::cout << "Test passed" << std::endl;

        // quad that is much bigger, intersection an ellipse
        Quad quad9(Point(-2.0*f, -2.0*f, 0.5*f + 0.3*2.0*f),
                   Point( 2.0*f, -2.0*f, 0.5*f - 0.3*2.0*f),
                   Point(-2.0*f,  2.0*f, 0.5*f + 0.3*2.0*f),
                   Point( 2.0*f,  2.0*f, 0.5*f - 0.3*2.0*f));
        result = intersect(cylSurface, quad9);
        if (result.size() > 1)
            throw std::runtime_error("9: Found more than one intersection");
        std::cout << "Checking test 9" << std::endl;
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipse); }, variant);
        if (!std::get<Ellipse>(result[0]).getPoint(0.0).isEqual(Point({-0.5*f, 0.0, 0.5*f + 0.3*0.5*f}), eps))
            throw std::runtime_error("Unexpected ellipse (test 9)");
        if (!std::get<Ellipse>(result[0]).getPoint(0.5).isEqual(Point({0.5*f, 0.0, 0.5*f - 0.3*0.5*f}), eps))
            throw std::runtime_error("Unexpected ellipse (test 9)");
        std::cout << "Test passed" << std::endl;

        // quad that touches in two points
        Quad quad10(Point(-0.5*f,  0.0*f, 0.5*f),
                    Point( 0.0*f, -0.25*f, 0.5*f),
                    Point( 0.0*f,  0.25*f, 0.5*f),
                    Point( 0.5*f,  0.0*f, 0.5*f));
        result = intersect(cylSurface, quad10);
        if (result.size() != 2)
            throw std::runtime_error("10: Did not find two intersections");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, variant);
        if (!std::get<Point>(result[0]).isEqual(Point({-0.5*f, 0.0, 0.5*f}), eps))
            throw std::runtime_error("Unexpected point 10");
        if (!std::get<Point>(result[1]).isEqual(Point({0.5*f, 0.0, 0.5*f}), eps))
            throw std::runtime_error("Unexpected point 10");
        std::cout << "Test passed" << std::endl;

        // quad that intersects in four points
        Quad quad11(Point(-0.5*f,  0.0*f, 0.5*f),
                    Point( 0.0*f, -0.5*f, 0.5*f),
                    Point( 0.0*f,  0.5*f, 0.5*f),
                    Point( 0.5*f,  0.0*f, 0.5*f));
        result = intersect(cylSurface, quad11);
        if (result.size() != 4)
            throw std::runtime_error("11: Did not find four intersections");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, variant);
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Point>(is).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected segment center (test 11)");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Point>(is).isEqual(Point(0.0*f, 0.5*f, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected segment center (test 11)");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Point>(is).isEqual(Point(-0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected segment center (test 11)");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Point>(is).isEqual(Point(0.0*f, -0.5*f, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected segment center (test 11)");
        std::cout << "Test passed" << std::endl;

        // quad that intersects in two short segments
        Quad quad12(Point(-0.6*f, 0.0*f, 0.5*f),
                    Point( 0.0*f, 0.0*f, 0.4*f),
                    Point( 0.0*f, 0.0*f, 0.6*f),
                    Point( 0.6*f, 0.0*f, 0.5*f));
        result = intersect(cylSurface, quad12);
        if (result.size() != 2)
            throw std::runtime_error("12: Did not find two intersections");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, variant);
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
        std::cout << "Test passed" << std::endl;

        // quad that intersects in two short arcs
        Quad quad13(Point(-0.6*f,  0.0*f, 0.5*f),
                    Point( 0.0*f, -0.1*f, 0.5*f),
                    Point( 0.0*f,  0.1*f, 0.5*f),
                    Point( 0.6*f,  0.0*f, 0.5*f));
        result = intersect(cylSurface, quad13);
        if (result.size() != 2)
            throw std::runtime_error("13: Did not find two intersections");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, variant);
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected ellipse arc (test 13)");
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error("Unexpected ellipse arc (test 13)");
        std::cout << "Test passed" << std::endl;

        // quad that intersects in two long arcs
        Quad quad14(Point(-0.6*f,  0.0*f, 0.5*f),
                    Point( 0.0*f, -0.4*f, 0.5*f),
                    Point( 0.0*f,  0.4*f, 0.5*f),
                    Point( 0.6*f,  0.0*f, 0.5*f));
        result = intersect(cylSurface, quad14);
        if (result.size() != 2)
            throw std::runtime_error("14: Did not find two intersections");
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, variant);
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
        std::cout << "Test passed" << std::endl;

        // quad that intersects in one segment and one point
        Quad quad15(Point(-0.6*f, 0.0*f, 0.5*f),
                    Point( 0.0*f, 0.0*f, 0.4*f),
                    Point( 0.0*f, 0.0*f, 0.6*f),
                    Point( 0.5*f, 0.0*f, 0.5*f));
        result = intersect(cylSurface, quad15);
        if (result.size() != 2)
            throw std::runtime_error("15: Did not find two intersections");
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result[0]);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result[1]);
        std::cout << "Test passed" << std::endl;

        // quad that intersects in one arc and one point
        Quad quad16(Point(-0.6*f,  0.0*f, 0.5*f),
                    Point( 0.0*f, -0.1*f, 0.5*f),
                    Point( 0.0*f,  0.1*f, 0.5*f),
                    Point( 0.5*f,  0.0*f, 0.5*f));
        result = intersect(cylSurface, quad16);
        if (result.size() != 2)
            throw std::runtime_error("16: Did not find two intersections");
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result[0]);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, result[1]);
        std::cout << "Test passed" << std::endl;

        // disk that intersects in one parabola
        Quad quad17(Point( 0.0*f, -1.0*f, 0.0*f),
                    Point( 2.0*f, -0.5*f, 0.5*f),
                    Point(-2.0*f, -0.5*f, 0.5*f),
                    Point( 0.0*f,  1.0*f, 2.0*f));
        result = intersect(cylSurface, quad17);
        if (result.size() != 1)
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
        std::cout << "Test passed" << std::endl;

        std::cout << "Test passed" << std::endl;

        std::cout << "All tests passed"  << std::endl;
    }

    return 0;
}
