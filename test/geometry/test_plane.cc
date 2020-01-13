#include <stdexcept>
#include <string>

#include <frackit/geometry/plane.hh>

//! tests some functionality of planes
int main()
{
    using ctype = double;

    using Plane3 = Frackit::Plane<ctype, 3>;
    using Direction3 = typename Plane3::Direction;
    using Point3 = typename Plane3::Point;

    Direction3 e1(Direction3::Vector({1.0, 0.0, 0.0}));
    Direction3 e2(Direction3::Vector({0.0, 1.0, 0.0}));
    Direction3 e3(Direction3::Vector({0.0, 0.0, 1.0}));

    // create three identical planes with different ctors
    Plane3 plane1(Point3(0.0, 0.0, 0.0), e1, e2, e3);
    Plane3 plane2(Point3(0.0, 0.0, 0.0), e3);
    Plane3 plane3(Point3(0.0, 0.0, 0.0), Point3(1.0, 0.0, 0.0), Point3(1.0, 1.0, 0.0));

    // project a segment on the planes
    typename Plane3::Segment s(Point3(0.0, 0.0, 0.0), Point3(1.0, 1.0, 1.0));
    const auto s1 = plane1.projection(s);
    const auto s2 = plane2.projection(s);
    const auto s3 = plane2.projection(s);

    // compare the resulting segments with the expected result
    Point3 source(0.0, 0.0, 0.0);
    Point3 target(1.0, 1.0, 0.0);

    using Vector3 = typename Direction3::Vector;
    if (Vector3(s1.source(), source).length() > 1e-7) throw std::runtime_error(std::string("Unexpected source 1"));
    if (Vector3(s1.target(), target).length() > 1e-7) throw std::runtime_error(std::string("Unexpected target 1"));

    if (Vector3(s2.source(), source).length() > 1e-7) throw std::runtime_error(std::string("Unexpected source 2"));
    if (Vector3(s2.target(), target).length() > 1e-7) throw std::runtime_error(std::string("Unexpected target 2"));

    if (Vector3(s3.source(), source).length() > 1e-7) throw std::runtime_error(std::string("Unexpected source 3"));
    if (Vector3(s3.target(), target).length() > 1e-7) throw std::runtime_error(std::string("Unexpected target 3"));

    const auto p1 = plane1.projection(target);
    const auto p2 = plane2.projection(target);
    const auto p3 = plane3.projection(target);

    if (Vector3(p1, Point3(1.0, 1.0, 0.0)).length() > 1e-7) throw std::runtime_error(std::string("Unexpected point projection 1"));
    if (Vector3(p2, Point3(1.0, 1.0, 0.0)).length() > 1e-7) throw std::runtime_error(std::string("Unexpected point projection 2"));
    if (Vector3(p3, Point3(1.0, 1.0, 0.0)).length() > 1e-7) throw std::runtime_error(std::string("Unexpected point projection 3"));

    std::cout << "All tests passed" << std::endl;
    return 0;
}
