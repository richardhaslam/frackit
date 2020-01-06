#include <frackit/geometry/disk.hh>
#include <frackit/geometry/ellipsearc.hh>
#include <frackit/geometry/cylindersurface.hh>

#include <frackit/magnitude/length.hh>
#include <frackit/intersection/intersect.hh>
#include <frackit/intersection/intersectionresult.hh>

enum IntersectionType { point, ellipseArc, ellipse, segment, empty };

template<class G>
void checkResultGeometry(const G& geometry, IntersectionType expected)
{
    throw std::runtime_error(std::string("Unexpected intersection geometry"));
}

void checkResultGeometry(const Frackit::EmptyIntersection<3>& empty, IntersectionType expected)
{
    std::cout << "Found empty intersection" << std::endl;
    if (expected != IntersectionType::empty)
        throw std::runtime_error(std::string("Unexpected empty intersection"));
}

template<class CT, int wd>
void checkResultGeometry(const Frackit::Point<CT, wd>& p, IntersectionType expected)
{
    std::cout << "Found intersection point at " << p << std::endl;
    if (expected != IntersectionType::point)
        throw std::runtime_error(std::string("Got an unexpected point intersection"));
}

template<class CT, int wd>
void checkResultGeometry(const Frackit::Segment<CT, wd>& segment, IntersectionType expected)
{
    std::cout << "Found intersection segment with corners "
              << segment.source() << " - " << segment.target()
              << " and length " << Frackit::computeLength(segment) << std::endl;
    if (expected != IntersectionType::segment)
        throw std::runtime_error(std::string("Got an unexpected segment intersection"));
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
        throw std::runtime_error(std::string("Got an unexpected arc intersection"));
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
        throw std::runtime_error(std::string("Got an unexpected ellipse intersection"));
}

template<class IS>
void checkResultGeometry(const std::vector<IS>& intersections, IntersectionType expected)
{
    for (const auto& is : intersections)
        checkResultGeometry(is, expected);
}

//! test disk-disk intersections
int main()
{
    using ctype = double;
    using CylinderSurface = Frackit::CylinderSurface<ctype>;
    using Segment = typename CylinderSurface::Segment;
    using Disk = Frackit::Disk<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Ellipse = typename Disk::Ellipse;
    using Vector = typename Direction::Vector;
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

        // disk that touches the cylinder surface in one point
        auto result = intersect(cylSurface, Disk(Point(-1.0*f, 0.0, 0.0), e1, e2, 0.5*f, 0.25*f));
        if (result.size() > 1)
            throw std::runtime_error(std::string("1: More than one intersection found"));
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, variant);
        const auto p1 = std::get<Point>(result[0]);
        if (!p1.isEqual(Point({-0.5*f, 0.0, 0.0}), eps))
            throw std::runtime_error(std::string("Unexpected point 1"));
        std::cout << "Test passed" << std::endl;

        // disk that describes a horizontal cross-section (circle) of the surface
        result = intersect(cylSurface, Disk(Point(0.0, 0.0, 0.5*f), e1, e2, 0.5*f, 0.5*f));
        if (result.size() > 1)
            throw std::runtime_error(std::string("2: More than one intersection found"));
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipse); }, variant);
        const auto ellipse1 = std::get<Ellipse>(result[0]);
        using std::abs;
        if (abs(Vector(ellipse1.normal())*e1) > eps)
            throw std::runtime_error(std::string("Unexpected normal direction of ellipse 1"));
        if (abs(Vector(ellipse1.normal())*e2) > eps)
            throw std::runtime_error(std::string("Unexpected normal direction of ellipse 1"));
        if (abs(1.0 - abs(Vector(ellipse1.normal())*e3)) > eps)
            throw std::runtime_error(std::string("Unexpected normal direction of ellipse 1"));
        std::cout << "Test passed" << std::endl;

        // disk that describes a vertical cross-section of the surface (two segments)
        result = intersect(cylSurface, Disk(Point(0.0, 0.0, 0.5*f), e1, e3, 2.0*f, 2.0*f));
        if (result.size() != 2)
            throw std::runtime_error(std::string("3: Did not find two intersection geometries"));
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, variant);
        if ( abs(Frackit::computeLength(std::get<Segment>(result[0])) - cylSurface.height()) > eps )
            throw std::runtime_error(std::string("Unexpected segment height (test 3)"));
        if ( abs(Frackit::computeLength(std::get<Segment>(result[1])) - cylSurface.height()) > eps )
            throw std::runtime_error(std::string("Unexpected segment height (test 3)"));
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Segment>(is).getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error(std::string("Unexpected segment (test 3)"));
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Segment>(is).getPoint(0.5).isEqual(Point(-0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error(std::string("Unexpected segment (test 3)"));
        std::cout << "Test passed" << std::endl;

        // disk that touches the cylinder surface in one point (opposite side)
        result = intersect(cylSurface, Disk(Point(1.0*f, 0.0, 0.0), e1, e2, 0.5*f, 0.25*f));
        if (result.size() > 1)
            throw std::runtime_error(std::string("4: More than one intersection found"));
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, variant);
        if (!std::get<Point>(result[0]).isEqual(Point({0.5*f, 0.0, 0.0}), eps))
            throw std::runtime_error(std::string("Unexpected point (test 4)"));
        std::cout << "Test passed" << std::endl;

        // no intersection
        result = intersect(cylSurface, Disk(Point(1.0*f + 1e-3*f, 0.0, 0.0), e1, e2, 0.5*f, 0.25*f));
        if (result.size() > 1)
            throw std::runtime_error(std::string("5: More than one intersection found"));
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::empty); }, variant);
        std::cout << "Test passed" << std::endl;

        // disk that intersects on one side only
        result = intersect(cylSurface, Disk(Point(1.0*f, 0.0, 0.5*f), e1, e2, 0.75*f, 0.25*f));
        if (result.size() > 1)
            throw std::runtime_error(std::string("6: More than one intersection found"));
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, variant);
        if (!std::get<EllipseArc>(result[0]).getPoint(0.5).isEqual(Point({0.5*f, 0.0, 0.5*f}), eps))
            throw std::runtime_error(std::string("Unexpected ellipse arc (test 6)"));
        std::cout << "Test passed" << std::endl;

        // disk that intersects on two sides
        result = intersect(cylSurface, Disk(Point(0.0, 0.0, 0.5*f), e1, e2, 0.75*f, 0.25*f));
        if (result.size() != 2)
            throw std::runtime_error(std::string("7: Did not find two intersections"));
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, variant);
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error(std::string("Unexpected ellipse arc (test 7)"));
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(-0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error(std::string("Unexpected ellipse arc (test 7)"));
        std::cout << "Test passed" << std::endl;

        // same as the previous test but with an inclined disk
        const auto e11 = Direction(Vector(1.0, 0.0, -0.3));
        const auto e21 = Direction(Vector(0.0, 1.0, 0.0));
        result = intersect(cylSurface, Disk(Point(0.0, 0.0, 0.5*f), e11, e21, 0.75*f, 0.25*f));
        if (result.size() != 2)
            throw std::runtime_error(std::string("8: Did not find two intersections"));
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, variant);
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f - 0.3*0.5*f), eps); }) != 1)
            throw std::runtime_error(std::string("Unexpected ellipse arc (test 8)"));
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(-0.5*f, 0.0, 0.5*f + 0.3*0.5*f), eps); }) != 1)
            throw std::runtime_error(std::string("Unexpected ellipse arc (test 8)"));
        std::cout << "Test passed" << std::endl;

        // disk that is much bigger, intersection an ellipse
        result = intersect(cylSurface, Disk(Point(0.0, 0.0, 0.5*f), e11, e21, 2.0*f, 2.0*f));
        if (result.size() > 1)
            throw std::runtime_error(std::string("9: Found more than one intersection"));
        std::cout << "Checking test 9" << std::endl;
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipse); }, variant);
        if (!std::get<Ellipse>(result[0]).getPoint(0.0).isEqual(Point({-0.5*f, 0.0, 0.5*f + 0.3*0.5*f}), eps))
            throw std::runtime_error(std::string("Unexpected ellipse (test 9)"));
        if (!std::get<Ellipse>(result[0]).getPoint(0.5).isEqual(Point({0.5*f, 0.0, 0.5*f - 0.3*0.5*f}), eps))
            throw std::runtime_error(std::string("Unexpected ellipse (test 9)"));
        std::cout << "Test passed" << std::endl;

        // disk that touches in two points
        result = intersect(cylSurface, Disk(Point(0.0, 0.0, 0.5*f), e1, e2, 0.5*f, 0.12*f));
        if (result.size() != 2)
            throw std::runtime_error(std::string("10: Did not find two intersections"));
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, variant);
        if (!std::get<Point>(result[0]).isEqual(Point({0.5*f, 0.0, 0.5*f}), eps))
            throw std::runtime_error(std::string("Unexpected point 10"));
        if (!std::get<Point>(result[1]).isEqual(Point({-0.5*f, 0.0, 0.5*f}), eps))
            throw std::runtime_error(std::string("Unexpected point 10"));
        std::cout << "Test passed" << std::endl;

        // disk that intersects in two short segments
        result = intersect(cylSurface, Disk(Point(0.0, 0.0, 0.5*f), e1, e3, 0.6*f, 0.12*f));
        if (result.size() != 2)
            throw std::runtime_error(std::string("11: Did not find two intersections"));
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, variant);
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Segment>(is).getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error(std::string("Unexpected segment center (test 11)"));
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<Segment>(is).getPoint(0.5).isEqual(Point(-0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error(std::string("Unexpected segment center (test 11)"));
        std::cout << "Test passed" << std::endl;

        // disk that intersects in two short arcs
        const auto e22 = Direction(Vector(0.0, 1.0, 0.3));
        result = intersect(cylSurface, Disk(Point(0.0, 0.0, 0.5*f), e1, e22, 0.6*f, 0.12*f));
        if (result.size() != 2)
            throw std::runtime_error(std::string("12: Did not find two intersections"));
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, variant);
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error(std::string("Unexpected ellipse arc (test 12)"));
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error(std::string("Unexpected ellipse arc (test 12)"));
        std::cout << "Test passed" << std::endl;

        // disk that intersects in two long arcs
        const auto e23 = Direction(Vector(0.0, 0.3, 1.0));
        result = intersect(cylSurface, Disk(Point(0.0, 0.0, 0.5*f), e1, e23, 2.0*f, 2.0*f));
        if (result.size() != 2)
            throw std::runtime_error(std::string("13: Did not find two intersections"));
        for (const auto& variant : result)
            std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, variant);
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error(std::string("Unexpected ellipse arc (test 13)"));
        if (std::count_if(result.begin(),
                          result.end(),
                          [f, eps] (const auto& is)
                          { return std::get<EllipseArc>(is).getPoint(0.5).isEqual(Point(0.5*f, 0.0, 0.5*f), eps); }) != 1)
            throw std::runtime_error(std::string("Unexpected ellipse arc (test 13)"));
        std::cout << "Test passed" << std::endl;

        // disk that intersects in one segment and one point
        result = intersect(cylSurface, Disk(Point(0.5*f, 0.0, 0.5*f), e1, e3, f, 0.2*f));
        if (result.size() != 2)
            throw std::runtime_error(std::string("14: Did not find two intersections"));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result[0]);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result[1]);
        std::cout << "Test passed" << std::endl;

        // disk that intersects in one arc and one point
        result = intersect(cylSurface, Disk(Point(0.5*f, 0.0, 0.5*f), e1, e22, f, 0.2*f));
        if (result.size() != 2)
            throw std::runtime_error(std::string("15: Did not find two intersections"));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result[0]);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, result[1]);
        std::cout << "Test passed" << std::endl;

        // disk that intersects in one parabola
        const auto e24 = Direction(Vector(0.0, 1.0, 0.75));
        result = intersect(cylSurface, Disk(Point(0.0, 0.0, 1.25*f), e1, e24, 2.0*f, 2.0*f));
        if (result.size() != 1)
            throw std::runtime_error(std::string("16: Did not find a single intersections"));
        if (abs(std::get<EllipseArc>(result[0]).source().z() - f) > eps)
            throw std::runtime_error(std::string("Unexpected ellipse arc"));
        if (abs(std::get<EllipseArc>(result[0]).target().z() - f) > eps)
            throw std::runtime_error(std::string("Unexpected ellipse arc"));
        if (abs(std::get<EllipseArc>(result[0]).getPoint(0.5).x()) > eps)
            throw std::runtime_error(std::string("Unexpected ellipse arc"));
        if (abs(std::get<EllipseArc>(result[0]).getPoint(0.5).y() + 0.5*f) > eps)
            throw std::runtime_error(std::string("Unexpected ellipse arc"));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, result[0]);
        std::cout << "Test passed" << std::endl;

        // disk that intersects in one parabola on its boundary
        const auto e25 = Direction(Vector(0.0, 1.0, 1.0));
        const auto e14 = Direction(Vector(-1.0, 0.0, 0.0));
        const auto cuttingDisk = Disk(Point(0.0, 0.0, 1.0*f), e25, e14, f*0.5/std::cos(M_PI/4.0), 0.5*f);
        result = intersect(cylSurface, cuttingDisk);
        if (result.size() != 1)
            throw std::runtime_error(std::string("17: Did not find a single intersections"));
        if (abs(std::get<EllipseArc>(result[0]).source().z() - f) > eps)
            throw std::runtime_error(std::string("Unexpected ellipse arc"));
        if (abs(std::get<EllipseArc>(result[0]).target().z() - f) > eps)
            throw std::runtime_error(std::string("Unexpected ellipse arc"));
        if (abs(std::get<EllipseArc>(result[0]).getPoint(0.5).x()) > eps)
            throw std::runtime_error(std::string("Unexpected ellipse arc"));
        if (abs(std::get<EllipseArc>(result[0]).getPoint(0.5).y() + 0.5*f) > eps)
            throw std::runtime_error(std::string("Unexpected ellipse arc"));
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::ellipseArc); }, result[0]);
        std::cout << "Test passed" << std::endl;

        // the arc is completely on the ellipse boundary
        const auto& arc = std::get<Frackit::EllipseArc<ctype, 3>>(result[0]);
        const auto& arcEllipse = arc.supportingEllipse();
        const auto& be = cuttingDisk.boundingEllipse();

        std::vector<ctype> params({0.0, 0.25, 0.5, 0.75, 1.0});
        for (auto param : params)
            if (!be.contains(arcEllipse.getPoint(param), eps))
                throw std::runtime_error(std::string("Point not on arc"));

        std::cout << "Test passed" << std::endl;

        std::cout << "All tests passed"  << std::endl;
    }

    return 0;
}
