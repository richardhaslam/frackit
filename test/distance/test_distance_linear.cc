#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

#include <frackit/geometry/point.hh>
#include <frackit/geometry/line.hh>
#include <frackit/geometry/vector.hh>
#include <frackit/distance/distance.hh>
#include <frackit/precision/precision.hh>

//! test some distance computation functionalities
int main()
{
    static constexpr int worldDim = 3;
    using ctype = double;
    using Point = Frackit::Point<ctype, worldDim>;
    using Line = Frackit::Line<ctype, worldDim>;
    using Vector = Frackit::Vector<ctype, worldDim>;
    using Segment = Frackit::Segment<ctype, worldDim>;

    using std::abs;
    using std::sqrt;
    std::vector<ctype> scales({1e-5, 1, 1e5});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;

        // distance between two points
        auto d = computeDistance(Point(0.0, 0.0, 0.0), Point(f, 0.0, 0.0));
        if ( abs(d - f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error("Point-Point distance wrong");

        // distance between line and point
        d = computeDistance(Point(0.0, 0.0, 0.0),
                            Line(Point(0.0, f, 0.0), Vector(0.0, 0.0, 1.0)));
        if ( abs(d - f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error("Point-Line distance wrong");

        // distance between segment and point (center)
        d = computeDistance(Point(0.0, 0.0, 0.0),
                            Segment(Point(0.0, -f, f), Point(0.0, f, f)));
        if ( abs(d - f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error("Point-Segment distance wrong (1)");

        // distance between segment and point (corner)
        d = computeDistance(Point(0.0, 0.0, 0.0),
                            Segment(Point(0.0, 0.0, f), Point(0.0, f, f)));
        if ( abs(d - f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error("Point-Segment distance wrong (2)");

        // distance between segment and point (corner)
        d = computeDistance(Point(0.0, 0.0, 0.0),
                            Segment(Point(0.0, f, f), Point(0.0, 2.0*f, f)));
        if ( abs(d - sqrt(2.0)*f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error("Point-Segment distance wrong (3)");
    }

    std::cout << "All tests passed" << std::endl;
    return 0;
}
