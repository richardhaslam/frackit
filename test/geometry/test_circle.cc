#include <stdexcept>
#include <string>

#include <frackit/geometry/circle.hh>

//! test some functionality of circles
int main()
{
    using ctype = double;

    using Circle = Frackit::Circle<ctype, 3>;
    using Point = typename Circle::Point;
    using Direction = typename Circle::Direction;
    using Vector = typename Direction::Vector;

    Circle circle(Point(0.0, 0.0, 0.0), Direction(Vector(0.0, 0.0, 1.0)), 1.0);

    // check normal vector
    using std::abs;
    if ( abs(abs(Vector(0.0, 0.0, 1.0)*Vector(circle.normal())) - 1.0) > 1e-7 )
        throw std::runtime_error(std::string("Unexpected circle normal"));

    // point that is just next to the circle
    if (circle.contains(Point(1.0 - 1e-6, 0.0, 0.0)))
        throw std::runtime_error(std::string("False positive contains() result"));

    // point that is just next to the circle
    if (circle.contains(Point(1.0 + 1e-6, 0.0, 0.0)))
        throw std::runtime_error(std::string("False positive contains() result"));

    // point that is on
    if (!circle.contains(Point(1.0, 0.0, 0.0)))
        throw std::runtime_error(std::string("False negative contains() result"));

    std::cout << "All tests passed" << std::endl;
    return 0;
}
