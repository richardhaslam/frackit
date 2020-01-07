#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

#include <frackit/geometry/disk.hh>
#include <frackit/distance/distance.hh>

//! test some functionality of ellipse arcs
int main()
{
    using ctype = double;
    using Disk = Frackit::Disk<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;

    // base directions
    Direction e1(Vector(1.0, 0.0, 0.0));
    Direction e2(Vector(0.0, 1.0, 0.0));
    Direction e3(Vector(0.0, 0.0, 1.0));

    using std::abs;
    using std::sqrt;
    std::vector<ctype> scales({1e-5, 1, 1e5});
    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;

        // distance between two parallel disks
        auto d = Frackit::computeDistance(Disk(Point(0.0, 0.0, 0.0), e1, e2, f, 0.5*f),
                                          Disk(Point(0.0, 0.0, f), e1, e2, f, 0.5*f));
        if ( abs(d - f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error(std::string("Test 1 failed"));

        d = Frackit::computeDistance(Disk(Point(0.0, 0.0, 0.0), e1, e2, f, 0.5*f),
                                     Disk(Point(f, 0.0, f), e1, e2, f, 0.5*f));
        if ( abs(d - f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error(std::string("Test 2 failed"));

        d = Frackit::computeDistance(Disk(Point(0.0, 0.0, 0.0), e1, e2, f, 0.5*f),
                                     Disk(Point(3.0*f, 0.0, f), e1, e2, f, 0.5*f));
        if ( abs(d - sqrt(2.0)*f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error(std::string("Test 3 failed"));

        // orthogonal disks being closest in the middle
        d = Frackit::computeDistance(Disk(Point(0.0, 0.0, 0.0), e1, e2, f, 0.5*f),
                                     Disk(Point(0.0, 0.0, 2.0*f), e1, e3, f, f));
        if ( abs(d - f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error(std::string("Test 4 failed"));

        d = Frackit::computeDistance(Disk(Point(0.0, 0.0, 0.0), e1, e2, f, 0.5*f),
                                     Disk(Point(0.0, 0.0, 1.0*f + 1e-3*f), e1, e3, f, f));
        if ( abs(d - 1e-3*f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error(std::string("Test 5 failed"));

        d = Frackit::computeDistance(Disk(Point(0.0, 0.0, 0.0), e1, e2, f, 0.5*f),
                                     Disk(Point(0.0, 0.0, 1.0*f), e1, e3, f, 0.5*f));
        if ( abs(d - 0.5*f) > Frackit::Precision<ctype>::confusion()*f )
            throw std::runtime_error(std::string("Test 6 failed"));
    }

    std::cout << "All tests passed" << std::endl;
    return 0;
}
