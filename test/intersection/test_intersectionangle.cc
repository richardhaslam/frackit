#include <cmath>
#include <variant>

#include <TopoDS_Edge.hxx>

#include <frackit/geometry/disk.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/cylinder.hh>
#include <frackit/geometryutilities/name.hh>
#include <frackit/intersection/intersect.hh>
#include <frackit/intersection/intersectionangle.hh>
#include <frackit/occ/breputilities.hh>
#include <frackit/common/math.hh>

// throw error that more than one intersection was found
template<class ISVariant>
void printIntersectionError(const std::vector<ISVariant>& is, const std::string& expected)
{
    if (is.size() == 0) throw std::runtime_error("No intersections found");

    std::cout << "The following intersections were found: "<< std::endl;
    auto printName = [] (const auto& g) { std::cout << Frackit::geometryName(g) << std::endl; };
    for (const auto& isGeom : is) std::visit(printName, isGeom);
    throw std::runtime_error(expected);
}

//! test intersection angles
int main()
{
    using namespace Frackit;
    using ctype = double;
    using Disk = Disk<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;
    using Quad = Quadrilateral<ctype, 3>;
    using Cylinder = Cylinder<ctype>;
    using Segment = typename Cylinder::Segment;

    std::vector<ctype> scales({1.0e-5, 1, 1e5});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;

        const Direction e1(Vector(1.0, 0.0, 0.0));
        const Direction e2(Vector(0.0, 1.0, 0.0));
        const Disk disk(Point(0.0, 0.0, 0.0), e1, e2, f, f);
        const auto diskShape = OCCUtilities::getShape(disk);

        using std::sin;
        using std::cos;
        using std::abs;

        static constexpr int numangles = 72; // 5° steps
        static constexpr ctype deltaAngle = 2.0*M_PI/numangles;

        // the computed angles are always in between 0 <= angle <= pi
        auto getExpectedAngle = [] (const auto angle)
        {
            using std::abs;
            if (angle > M_PI/2.0 - Frackit::Precision<ctype>::angular()
                && angle < M_PI + Frackit::Precision<ctype>::angular())
                return abs(M_PI - angle);

            if (angle > M_PI - Frackit::Precision<ctype>::angular()
                && angle < 3.0/2.0*M_PI + Frackit::Precision<ctype>::angular())
                return angle - M_PI;

            if (angle > 3.0/2.0*M_PI - Frackit::Precision<ctype>::angular()
                && angle < 2.0*M_PI + Frackit::Precision<ctype>::angular())
                return abs(2.0*M_PI - angle);

            return angle;
        };

        std::cout << " -- testing segment-like disk-disk and face-face intersections" << std::endl;;
        for (int i = 0; i < numangles; ++i)
        {
            ctype angle = i*deltaAngle;
            const Direction majAxis( Vector(cos(angle), 0.0, sin(angle)) );
            const Disk tmp(Point(0.0, 0.0, 0.0), majAxis, e2, 0.5*f, 0.5*f);
            const auto tmpShape = OCCUtilities::getShape(tmp);
            const auto expectedAngle = getExpectedAngle(angle);

            // skip face intersections
            if ( abs(angle) < Frackit::Precision<ctype>::angular()
                 || abs(angle - M_PI) < Frackit::Precision<ctype>::angular() )
                continue;

            IntersectionAngle<ctype> angleEngine;
            const auto isShapes = intersect(diskShape, tmpShape);
            if (isShapes.size() != 1) printIntersectionError(isShapes, "Expected segment intersection (from shapes)");
            if (!std::holds_alternative<TopoDS_Edge>(isShapes[0])) printIntersectionError(isShapes, "Expected segment intersection (from shapes)");

            const auto compAngleShapes = angleEngine(diskShape, tmpShape, isShapes);
            if (abs(compAngleShapes - expectedAngle) > Frackit::Precision<ctype>::angular())
            {
                std::cout << "computed vs given angle (in degrees): "
                          << toDegrees(compAngleShapes) << " - " << toDegrees(expectedAngle) << std::endl;
                throw std::runtime_error("Wrong computed angle with disk shapes");
            }

            const auto isGeom = intersect(disk, tmp);
            if (!std::holds_alternative<Segment>(isGeom)) printIntersectionError(isShapes, "Expected segment intersection (from internal geoms)");
            const auto compAngleDisks = angleEngine(disk, tmp, isGeom);
            if (abs(compAngleDisks - expectedAngle) > Frackit::Precision<ctype>::angular())
            {
                std::cout << "computed vs given angle (in degrees): "
                          << toDegrees(compAngleDisks) << " - " << toDegrees(expectedAngle) << std::endl;
                throw std::runtime_error("Wrong computed angle with internal geometries");
            }
        }

        // test angles at point intersections for disk-quad & face-face
        std::cout << " -- testing point-like disk-quad and face-face intersections" << std::endl;;
        for (int i = 0; i < numangles; ++i)
        {
            ctype angle = i*deltaAngle;
            const Quad quad(Point(0.0, 0.0, 0.0),
                            Point(f, cos(angle)*0.8*f, sin(angle)*0.8*f),
                            Point(-f, cos(angle)*f, sin(angle)*f),
                            Point( f, cos(angle)*f, sin(angle)*f));

            const auto quadShape = OCCUtilities::getShape(quad);
            const auto expectedAngle = getExpectedAngle(angle);

            // skip zero angles (in which case they would intersect in a face)
            if ( abs(angle) < Frackit::Precision<ctype>::angular()
                 || abs(angle - M_PI) < Frackit::Precision<ctype>::angular() )
                continue;

            IntersectionAngle<ctype> angleEngine;
            const auto isShapes = intersect(diskShape, quadShape);
            if (isShapes.size() != 1) printIntersectionError(isShapes, "Expected point intersection (from shapes)");
            if (!std::holds_alternative<Point>(isShapes[0])) printIntersectionError(isShapes, "Expected point intersection (from shapes)");

            const auto compAngleShapes = angleEngine(diskShape, quadShape, isShapes);
            if (abs(compAngleShapes - expectedAngle) > Frackit::Precision<ctype>::angular())
            {
                std::cout << "computed vs given angle (in degrees): "
                          << toDegrees(compAngleShapes) << " - " << toDegrees(expectedAngle) << std::endl;
                throw std::runtime_error("Wrong computed angle with disk/quad shapes");
            }

            const auto isGeom = intersect(disk, quad);
            if (!std::holds_alternative<Point>(isGeom)) printIntersectionError(isShapes, "Expected point intersection (from internal geoms)");

            const auto compAngleGeoms = angleEngine(disk, quad, isGeom);
            if (abs(compAngleGeoms - expectedAngle) > Frackit::Precision<ctype>::angular())
            {
                std::cout << "computed vs given angle (in degrees): "
                          << toDegrees(compAngleGeoms) << " - " << toDegrees(expectedAngle) << std::endl;
                throw std::runtime_error("Wrong computed angle with internal geometries");
            }
        }

        // test cylinder surface - planar geometry intersection angle
        // test only between -pi/2 and pi/2, as otherwise we wouldn't have point intersection
        std::cout << " -- testing point-like cylinder surface with quad/face intersections" << std::endl;
        for (int i = 0; i < numangles/2; ++i)
        {
            ctype angle = i*deltaAngle;
            const auto expectedAngle = getExpectedAngle(angle);

            // construct cylinder and quad such that it touches the cylinder at the desired angle
            const Cylinder cyl(f, f, -0.5*f);
            const Quad quad(Point(f, -f, 0.0), Point(f,  f, 0.0),
                            Point(f + abs(sin(angle)*f), -f, cos(angle)*f),
                            Point(f + abs(sin(angle)*f),  f, cos(angle)*f));

            const auto cylSurface = cyl.lateralFace();
            const auto cylSurfaceShape = OCCUtilities::getShape(cylSurface);
            const auto quadShape = OCCUtilities::getShape(quad);

            // skip zero angles (in which case they would intersect in a segment, is tested afterwards)
            if ( abs(angle) < Frackit::Precision<ctype>::angular()
                 || abs(abs(angle) - M_PI) < Frackit::Precision<ctype>::angular() )
                continue;

            IntersectionAngle<ctype> angleEngine;
            const auto isShapes = intersect(quadShape, cylSurfaceShape);
            if (isShapes.size() != 1) printIntersectionError(isShapes, "Expected point intersection (from shapes)");
            if (!std::holds_alternative<Point>(isShapes[0])) printIntersectionError(isShapes, "Expected point intersection (from shapes)");

            const auto compAngleShapes = angleEngine(quadShape, cylSurfaceShape, isShapes);
            if (abs(compAngleShapes - expectedAngle) > Frackit::Precision<ctype>::angular())
            {
                std::cout << "computed vs given angle (in degrees): "
                          << toDegrees(compAngleShapes) << " - " << toDegrees(expectedAngle) << std::endl;
                throw std::runtime_error("Wrong computed angle with quad/cylinder surface shapes");
            }

            const auto isGeoms = intersect(quad, cylSurface);
            if (isGeoms.size() != 1) printIntersectionError(isGeoms, "Expected point intersection (from internal geoms)");
            if (!std::holds_alternative<Point>(isGeoms[0])) printIntersectionError(isGeoms, "Expected point intersection (from internal geoms)");

            const auto compAngleGeoms = angleEngine(quad, cylSurface, isGeoms);
            if (abs(compAngleGeoms - expectedAngle) > Frackit::Precision<ctype>::angular())
            {
                std::cout << "computed vs given angle (in degrees): "
                          << toDegrees(compAngleGeoms) << " - " << toDegrees(expectedAngle) << std::endl;
                throw std::runtime_error("Wrong computed angle at cylinder with internal geometries");
            }
        }

        // test cylinder surface - planar geometry intersection angle
        // test only between -pi/2 and pi/2, as otherwise we maybe wouldn't have segment intersection
        std::cout << " -- testing segment-like cylinder surface with quad/face intersections" << std::endl;
        for (int i = 0; i < numangles/2; ++i)
        {
            ctype angle = i*deltaAngle;
            const auto expectedAngle = getExpectedAngle(angle);

            // construct cylinder and quad such that it touches the cylinder at the desired angle
            const Cylinder cyl(f, f, -0.5*f);
            const Quad quad(Point(0.0,          f,                     -f/4.0),
                            Point(0.0,          f,                      f/4.0),
                            Point(cos(angle)*f, f + abs(sin(angle))*f, -f/4.0),
                            Point(cos(angle)*f, f + abs(sin(angle))*f,  f/4.0));

            const auto cylSurface = cyl.lateralFace();
            const auto cylSurfaceShape = OCCUtilities::getShape(cylSurface);
            const auto quadShape = OCCUtilities::getShape(quad);

            // skip zero angles (in which case they would intersect in a segment, is tested afterwards)
            if ( abs(angle) < Frackit::Precision<ctype>::angular()
                 || abs(abs(angle) - M_PI) < Frackit::Precision<ctype>::angular() )
                continue;

            IntersectionAngle<ctype> angleEngine;
            const auto isShapes = intersect(quadShape, cylSurfaceShape);
            if (isShapes.size() != 1) printIntersectionError(isShapes, "Expected segment intersection (from shapes)");
            if (!std::holds_alternative<TopoDS_Edge>(isShapes[0])) printIntersectionError(isShapes, "Expected segment intersection (from shapes)");

            const auto compAngleShapes = angleEngine(quadShape, cylSurfaceShape, isShapes);
            if (abs(compAngleShapes - expectedAngle) > Frackit::Precision<ctype>::angular())
            {
                std::cout << "computed vs given angle (in degrees): "
                          << toDegrees(compAngleShapes) << " - " << toDegrees(expectedAngle) << std::endl;
                throw std::runtime_error("Wrong computed angle with quad/cylinder surface shapes");
            }

            const auto isGeoms = intersect(quad, cylSurface);
            if (isGeoms.size() != 1) printIntersectionError(isGeoms, "Expected segment intersection (from internal geoms)");
            if (!std::holds_alternative<Segment>(isGeoms[0])) printIntersectionError(isGeoms, "Expected segment intersection (from internal geoms)");

            const auto compAngleGeoms = angleEngine(quad, cylSurface, isGeoms);
            if (abs(compAngleGeoms - expectedAngle) > Frackit::Precision<ctype>::angular())
            {
                std::cout << "computed vs given angle (in degrees): "
                          << toDegrees(compAngleGeoms) << " - " << toDegrees(expectedAngle) << std::endl;
                throw std::runtime_error("Wrong computed angle at cylinder with internal geometries");
            }
        }

        std::cout << "All tests passed"  << std::endl;
    }

    return 0;
}
