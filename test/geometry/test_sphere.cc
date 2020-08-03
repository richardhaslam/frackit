#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

#include <frackit/geometry/sphere.hh>
#include <frackit/precision/defaultepsilon.hh>

//! test some functionality of triangles
int main()
{
    using ctype = double;
    using Sphere = Frackit::Sphere<ctype>;
    using Point = typename Sphere::Point;

    std::vector<ctype> scales({1e-5, 1, 1e5});
    for (auto f : scales)
    {
        Sphere sphere1(Point(0.0, 0.0, 0.0), f);
        Sphere sphere2(Point(2.0*f, 0.0, 0.0), f);
        const auto eps = defaultEpsilon(sphere1);

        using std::abs;
        if ( abs(sphere1.volume() - sphere2.volume()) > eps*eps*eps )
            throw std::runtime_error("Volumes are not the same");

        if ( abs(sphere1.volume() - 4.0/3.0*M_PI*f*f*f) > eps*eps*eps)
            throw std::runtime_error("Unexpected sphere volume");

        // check contains() query
        if (!sphere1.contains(Point(f, 0.0, 0.0), eps))
            throw std::runtime_error("contains() query 1 failed");
        if (!sphere1.contains(Point(0.0, 0.5*f - 1e-6*f, 0.0), eps))
            throw std::runtime_error("contains() query 2 failed");
        if (sphere2.contains(Point(f - 1e-6*f, 0.0, 0.0), eps))
            throw std::runtime_error("contains() query 3 failed");
        if (!sphere1.contains(Point(0, 0, f), eps))
            throw std::runtime_error("contains() query 4 failed");
        if (!sphere2.contains(Point(2.0*f, 0, f-1e-6*f), eps))
            throw std::runtime_error("contains() query 5 failed");
    }

    std::cout << "All tests passed" << std::endl;
    return 0;
}
