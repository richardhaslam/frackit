#include <frackit/geometry/point.hh>
#include <frackit/geometry/vector.hh>

//! tests some functionality of points
int main()
{
    using ctype = double;

    Frackit::Point<ctype, 1> p1(1.0);
    Frackit::Point<ctype, 2> p2(1.0, 2.0);
    Frackit::Point<ctype, 3> p3(1.0, 2.0, 3.0);

    std::cout << "P1: " << p1 << std::endl;
    std::cout << "P2: " << p2 << std::endl;
    std::cout << "P3: " << p3 << std::endl;

    // test 3d vector addition
    const auto p32 = p3 + Frackit::Vector<ctype, 3>(1.0, 2.0, 3.0);
    p3 += Frackit::Vector<ctype, 3>(1.0, 2.0, 3.0);
    if (!p3.isEqual(p32))
        throw std::runtime_error(std::string("p32 != p3"));

    return 0;
}
