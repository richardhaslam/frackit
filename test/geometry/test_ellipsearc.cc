#include <stdexcept>
#include <string>
#include <cmath>
#include <vector>

#include <frackit/magnitude/length.hh>
#include <frackit/precision/precision.hh>
#include <frackit/geometry/ellipsearc.hh>

//! test some functionality of ellipse arcs
int main()
{
    using ctype = double;

    using EllipseArc = Frackit::EllipseArc<ctype, 3>;
    using Ellipse = typename EllipseArc::Ellipse;
    using Point = typename EllipseArc::Point;
    using Direction = typename EllipseArc::Direction;
    using Vector = typename Direction::Vector;

    std::vector<ctype> scales({1e-5, 1, 1e5});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;
        const auto majAxis = f;
        const auto minAxis = 0.5*f;

        Vector e1(1.0, 0.0, 0.0);
        Vector e2(0.0, 1.0, 0.0);
        Ellipse ellipse(Point(0.0, 0.0, 0.0), e1, e2, majAxis, minAxis);

        const auto p = ellipse.getPointFromAngle(M_PI/4.0);
        const auto source = ellipse.getPointFromAngle(0.0);
        const auto target = ellipse.getPointFromAngle(M_PI/2.0);
        EllipseArc ellipseArc(ellipse, source, target);

        // check validity of source and target points
        if (!source.isEqual(Point(majAxis, 0.0, 0.0), f*Frackit::Precision<ctype>::confusion()))
            throw std::runtime_error("Unexpected source point");
        if (!target.isEqual(Point(0.0, minAxis, 0.0), f*Frackit::Precision<ctype>::confusion()))
            throw std::runtime_error("Unexpected source point");

        using std::abs;
        if ( abs(ellipseArc.getAngle(source)) > Frackit::Precision<ctype>::angular() )
            throw std::runtime_error("Wrong source angle");
        if ( abs(ellipseArc.getAngle(target) - M_PI/2.0) > Frackit::Precision<ctype>::angular() )
            throw std::runtime_error("Wrong target angle");
        if ( abs(ellipseArc.getAngle(p) - M_PI/4.0) > Frackit::Precision<ctype>::angular() )
            throw std::runtime_error("Wrong point angle");
        if ( abs(ellipseArc.getParam(p) - 1.0/8.0) > Frackit::Precision<ctype>::confusion() )
            throw std::runtime_error("Wrong point param");

        // check contains() queries
        if (!ellipseArc.contains(source))
            throw std::runtime_error("Source contains() query failed");
        if (!ellipseArc.contains(source, false))
            throw std::runtime_error("Source contains() query failed");

        if (!ellipseArc.contains(target))
            throw std::runtime_error("Target contains() query failed");
        if (!ellipseArc.contains(target, false))
            throw std::runtime_error("Target contains() query failed");

        if (!ellipseArc.contains(p))
            throw std::runtime_error("Point contains() query failed");
        if (!ellipseArc.contains(p, false))
            throw std::runtime_error("Point contains() query failed");

        // check points which are slightly outside
        const auto p2 = ellipse.getPointFromAngle(M_PI/2.0 + 1e-4);
        const auto p3 = ellipse.getPointFromAngle(2.0*M_PI - 1e-4);
        if (ellipseArc.contains(p2))
            throw std::runtime_error("P2 contains() query failed");
        if (ellipseArc.contains(p2, false))
            throw std::runtime_error("P2 contains() query failed");
        if (ellipseArc.contains(p3))
            throw std::runtime_error("P3 contains() query failed");
        if (ellipseArc.contains(p3, false))
            throw std::runtime_error("P3 contains() query failed");

        // check points which are slightly outside
        if (!ellipseArc.contains(ellipseArc.getPoint(0.0)))
            throw std::runtime_error("contains() query with computed point 1 failed");
        if (!ellipseArc.contains(ellipseArc.getPoint(0.25)))
            throw std::runtime_error("contains() query with computed point 2 failed");
        if (!ellipseArc.contains(ellipseArc.getPoint(0.5)))
            throw std::runtime_error("contains() query with computed point 3 failed");
        if (!ellipseArc.contains(ellipseArc.getPoint(0.75)))
            throw std::runtime_error("contains() query with computed point 4 failed");
        if (!ellipseArc.contains(ellipseArc.getPoint(1.0)))
            throw std::runtime_error("contains() query with computed point 5 failed");

        // create arcs describing a full ellipse
        EllipseArc fullArc1(ellipse, source, source);
        if (!fullArc1.isFullEllipse())
            throw std::runtime_error("Arc1 is not a full ellipse");
        EllipseArc fullArc2(ellipse, p, p);
        if (!fullArc2.isFullEllipse())
            throw std::runtime_error("Arc2 is not a full ellipse");

        // check arc length computation
        const ctype eps = f*Frackit::Precision<ctype>::confusion();
        if ( abs(Frackit::computeLength(fullArc1) - Frackit::computeLength(ellipse)) > eps )
            throw std::runtime_error("Full arc 1 does not have correct length");
        if ( abs(Frackit::computeLength(fullArc2) - Frackit::computeLength(ellipse)) > eps )
            throw std::runtime_error("Full arc 2 does not have correct length");

        // create circular arcs to test lengths on various orientations
        std::vector<Ellipse> circles;
        circles.emplace_back(Point(0.0, 0.0, 0.0), e1, e2, f, f);
        circles.emplace_back(Point(0.0, 0.0, 0.0), e1, Vector(0.0, 1.0, 1.0), f, f);
        circles.emplace_back(Point(0.0, 0.0, 0.0), Vector(1.0, 0.0, 1.0), e2, f, f);
        circles.emplace_back(Point(0.0, 0.0, 0.0), Vector(1.0, 1.0, 1.0), Vector(-1.0, 1.0, 0.0), f, f);

        std::vector<ctype> angles({0.0, M_PI/4.0, M_PI/2.0, 3.0*M_PI/4.0, M_PI, 3.0*M_PI/2.0});
        for (const auto& c : circles)
        {
            for (auto a1 : angles)
                for (auto a2 : angles)
                {
                    if (a1 == a2) continue;

                    const auto arc = EllipseArc(c, c.getPointFromAngle(a1), c.getPointFromAngle(a2));
                    const auto deltaAngle = a2 < a1 ? 2.0*M_PI - (a1 - a2) : a2 - a1;;

                    if ( abs(Frackit::computeLength(arc) - deltaAngle*f) > eps )
                        throw std::runtime_error("Wrong arc length");
                }
        }
    }

    std::cout << "All tests passed" << std::endl;
    return 0;
}
