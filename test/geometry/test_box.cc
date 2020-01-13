#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

#include <frackit/geometry/box.hh>
#include <frackit/precision/defaultepsilon.hh>

//! test some functionality of triangles
int main()
{
    using ctype = double;
    using Box = Frackit::Box<ctype>;
    using Point = typename Box::Point;

    std::vector<ctype> scales({1e-5, 1, 1e5});
    for (auto f : scales)
    {
        Box box(0.0, 0.0, 0.0, f, f, f);
        const auto eps = defaultEpsilon(box);

        using std::abs;
        if ( abs(box.volume() - f*f*f) > eps*eps*eps)
            throw std::runtime_error(std::string("Unexpected box volume"));

        // check edge lengths
        for (unsigned int i = 0; i < box.numEdges(); ++i)
            if ( abs(box.edge(i).length() - f) > eps)
                throw std::runtime_error(std::string("Edge " + std::to_string(i) + " has unexpected length"));

        // check face areas
        for (unsigned int i = 0; i < box.numFaces(); ++i)
            if ( abs(box.face(i).area() - f*f) > eps*eps)
                throw std::runtime_error(std::string("Face " + std::to_string(i) + " has unexpected area"));

        // check contains() query
        if (!box.contains(Point(0.5*f, 0.0, 0.0), eps))
            throw std::runtime_error(std::string("contains() query 1 failed"));
        if (!box.contains(Point(0.5*f, 1e-6*f, 0.0), eps))
            throw std::runtime_error(std::string("contains() query 2 failed"));
        if (box.contains(Point(0.5*f, -1e-6*f, 0.0), eps))
            throw std::runtime_error(std::string("contains() query 3 failed"));
        if (!box.contains(Point(0.5*f, 0.5*f, f), eps))
            throw std::runtime_error(std::string("contains() query 1 failed"));
        if (!box.contains(Point(0.5*f, 0.5*f, f-1e-6*f), eps))
            throw std::runtime_error(std::string("contains() query 2 failed"));
        if (box.contains(Point(0.5*f, 0.5*f, f+1e-6*f), eps))
            throw std::runtime_error(std::string("contains() query 3 failed"));
    }

    std::cout << "All tests passed" << std::endl;
    return 0;
}