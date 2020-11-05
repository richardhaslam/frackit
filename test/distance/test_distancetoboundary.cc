#include <cmath>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include <frackit/geometry/disk.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/polygon.hh>

#include <frackit/distance/distancetoboundary.hh>
#include <frackit/precision/precision.hh>
#include <frackit/intersection/intersect.hh>

//! test distance to boundary functionality
int main()
{
    using ctype = double;
    using Disk = Frackit::Disk<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;
    using Segment = Frackit::Segment<ctype, 3>;
    using Quad = Frackit::Quadrilateral<ctype, 3>;
    using Polygon = Frackit::Polygon<ctype, 3>;

    // base directions
    Direction e1(Vector(1.0, 0.0, 0.0));
    Direction e2(Vector(0.0, 1.0, 0.0));
    Direction e3(Vector(0.0, 0.0, 1.0));

    // lambda to print distance and throw error
    auto printAndThrow = [] (const auto d, const std::string& msg)
    {
        std::cout << "Computed a distance of " << d << std::endl;
        throw std::runtime_error(msg);
    };

    using std::abs;
    std::vector<ctype> scales({1e-3, 1, 1e3});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;
        const auto disk = Disk(Point(0.0, 0.0, 0.0), e1, e2, f, 0.5*f);
        const auto eps = Frackit::Precision<ctype>::confusion()*f;

        auto d = computeDistanceToBoundary(Point(0.0, 0.0, 0.0), disk);
        if (abs(d - 0.5*f) > eps) printAndThrow(d, "Test 1 failed");

        d = computeDistanceToBoundary(Point(f-1e-6*f, 0.0, 0.0), disk);
        if (abs(d - 1e-6*f) > eps) printAndThrow(d, "Test 2 failed");

        d = computeDistanceToBoundary(Point(f, 0.0, 1e-6*f), disk);
        if (abs(d - 1e-6*f) > eps) printAndThrow(d, "Test 3 failed");

        const Segment s1(Point(0.0, 0.0, 0.0), Point(f-1e-6*f, 0.0, 0.0));
        d = computeDistanceToBoundary(s1, disk);
        if (abs(d - 1e-6*f) > eps) printAndThrow(d, "Test 4 failed");

        // intersect with other disk (gives TopoDS_Edge) and try again
        const auto disk2 = Disk(Point(0.0, 0.0, 0.0), e2, e3, 0.5*f-1e-6*f, 0.25*f);
        const auto is = intersect(disk, disk2, eps);
        if (!std::holds_alternative<Segment>(is)) throw std::runtime_error("Expected segment intersection");
        d = computeDistanceToBoundary(std::get<Segment>(is), disk);
        if (abs(d - 1e-6*f) > eps) printAndThrow(d, "Test 5 failed");

        // do the same with a quadrilateral
        const Quad quad(Point(-0.5*f, -0.5*f, 0.0), Point(0.5*f, -0.5*f, 0.0),
                        Point(-0.5*f,  0.5*f, 0.0), Point(0.5*f,  0.5*f, 0.0));
        const auto is2 = intersect(quad, disk2, eps);
        if (!std::holds_alternative<Segment>(is2)) throw std::runtime_error("Expected segment intersection");
        d = computeDistanceToBoundary(Frackit::OCCUtilities::getShape(std::get<Segment>(is2)), quad);
        if (abs(d - 1e-6*f) > eps) printAndThrow(d, "Test 6 failed");

        // do the same with the quadrilateral represented by a polygon
        using CornerVector = std::vector<Point>;
        const Polygon poly(CornerVector{{Point(-0.5*f, -0.5*f, 0.0), Point(0.5*f, -0.5*f, 0.0),
                                         Point(0.5*f,  0.5*f, 0.0), Point(-0.5*f,  0.5*f, 0.0)}});
        const auto is3 = intersect(poly, disk2, eps);
        if (!std::holds_alternative<Segment>(is3)) throw std::runtime_error("Expected segment intersection");
        d = computeDistanceToBoundary(Frackit::OCCUtilities::getShape(std::get<Segment>(is3)), poly);
        if (abs(d - 1e-6*f) > eps) printAndThrow(d, "Test 7 failed");
    }

    std::cout << "All tests passed" << std::endl;
    return 0;
}
