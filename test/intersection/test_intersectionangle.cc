#include <cmath>

#include <frackit/geometry/disk.hh>
#include <frackit/intersection/intersect.hh>
#include <frackit/intersection/intersectionangle.hh>
#include <frackit/occ/breputilities.hh>
#include <frackit/common/math.hh>

//! test intersection angles
int main()
{
    using namespace Frackit;
    using ctype = double;
    using Disk = Disk<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;

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

        // test angles at segment intersections
        static constexpr int numangles = 72; // 5° steps
        static constexpr ctype deltaAngle = 2.0*M_PI/numangles;
        for (int i = 0; i < numangles; ++i)
        {
            ctype angle = i*deltaAngle;
            const Direction majAxis( Vector(cos(angle), 0.0, sin(angle)) );
            const Disk tmp(Point(0.0, 0.0, 0.0), majAxis, e2, 0.5*f, 0.5*f);
            const auto tmpShape = OCCUtilities::getShape(tmp);

            // computed angles are below 90°
            if (angle > M_PI/2.0 - Frackit::Precision<ctype>::angular()
                && angle < M_PI + Frackit::Precision<ctype>::angular())
                angle = abs(M_PI - angle);

            if (angle > M_PI - Frackit::Precision<ctype>::angular()
                && angle < 3.0/2.0*M_PI + Frackit::Precision<ctype>::angular())
                angle -= M_PI;

            if (angle > 3.0/2.0*M_PI - Frackit::Precision<ctype>::angular()
                && angle < 2.0*M_PI + Frackit::Precision<ctype>::angular())
                angle = abs(2.0*M_PI - angle);

            IntersectionAngle<ctype> angleEngine;
            const auto compAngleShapes = angleEngine(diskShape, tmpShape, intersect(diskShape, tmpShape));
            if (abs(compAngleShapes - angle) > Frackit::Precision<ctype>::angular())
            {
                std::cout << "computed vs given angle (in degrees): "
                          << toDegrees(compAngleShapes) << " - " << toDegrees(angle) << std::endl;
                throw std::runtime_error("Wrong computed angle with disk shapes");
            }

            // skip face intersections (currently not implenented)
            if ( angle < Frackit::Precision<ctype>::angular()
                 || abs(angle - M_PI) < Frackit::Precision<ctype>::angular() )
                continue;

            const auto compAngleDisks = angleEngine(disk, tmp, intersect(disk, tmp));
            if (abs(compAngleDisks - angle) > Frackit::Precision<ctype>::angular())
            {
                std::cout << "computed vs given angle (in degrees): "
                          << toDegrees(compAngleDisks) << " - " << toDegrees(angle) << std::endl;
                throw std::runtime_error("Wrong computed angle with disk shapes");
            }
        }

        std::cout << "All tests passed"  << std::endl;
    }

    return 0;
}
