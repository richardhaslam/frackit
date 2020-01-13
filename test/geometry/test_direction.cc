#include <stdexcept>
#include <string>

#include <frackit/geometry/direction.hh>

//! tests some functionality of directions
int main()
{
    using ctype = double;

    using Direction1 = Frackit::Direction<ctype, 1>;
    using Direction2 = Frackit::Direction<ctype, 2>;
    using Direction3 = Frackit::Direction<ctype, 3>;

    using Vector1 = typename Direction1::Vector;
    using Vector2 = typename Direction2::Vector;
    using Vector3 = typename Direction3::Vector;

    Vector1 v1(Vector1::Point(1.0), Vector1::Point(3.0));
    Vector2 v2(Vector2::Point(1.0, 2.0), Vector2::Point(2.0, 4.0));
    Vector3 v3(Vector3::Point(1.0, 2.0, 3.0), Vector3::Point(2.0, 4.0, 6.0));

    Direction1 d1(v1);
    Direction2 d2(v2);
    Direction3 d3(v3);

    using std::abs;
    using std::pow;
    using std::sqrt;

    const auto norm1 = abs(d1.x());
    const auto norm2 = sqrt(d2.x()*d2.x()+d2.y()*d2.y());
    const auto norm3 = sqrt(d3.x()*d3.x()+d3.y()*d3.y()+d3.z()*d3.z());

    if (abs(1.0 - norm1) > 1e-6) throw std::runtime_error(std::string("1d direction is " + std::to_string(norm1) + ", not unity"));
    if (abs(1.0 - norm2) > 1e-6) throw std::runtime_error(std::string("2d direction is " + std::to_string(norm2) + ", not unity"));
    if (abs(1.0 - norm3) > 1e-6) throw std::runtime_error(std::string("3d direction is " + std::to_string(norm3) + ", not unity"));

    std::cout << d1 << std::endl;
    std::cout << d2 << std::endl;
    std::cout << d3 << std::endl;

    return 0;
}
