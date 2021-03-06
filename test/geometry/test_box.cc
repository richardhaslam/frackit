#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

#include <frackit/geometry/point.hh>
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/vector.hh>
#include <frackit/geometry/box.hh>
#include <frackit/precision/defaultepsilon.hh>

//! test some functionality of triangles
int main()
{
    using ctype = double;
    using Box = Frackit::Box<ctype>;
    using Point = typename Box::Point;
    using Segment = typename Box::Segment;
    using Vector = Frackit::Vector<ctype, 3>;

    std::vector<ctype> scales({1e-5, 1, 1e5});
    for (auto f : scales)
    {
        Box box(0.0, 0.0, 0.0, f, f, f);
        const auto eps = defaultEpsilon(box);

        // also, test construction from corner points
        Point pMin(0.0, 0.0, 0.0);
        Point pMax(f,   f,   f  );
        Box box2(pMin, pMax);

        using std::abs;
        for (unsigned int i = 0; i < box.numCorners(); ++i)
            if (!box.corner(i).isEqual(box2.corner(i), eps))
                throw std::runtime_error("Construction from point led to different box");

        const auto diag = Vector(box.corner(0), box.corner(7)).length();
        const auto diag2 = Vector(box2.corner(0), box2.corner(7)).length();
        if ( abs(diag - diag2) > eps ) throw std::runtime_error("Box diagonals not equal");
        if ( abs(diag - sqrt(3.0)*f) > eps ) throw std::runtime_error("Unexpected box diagonal");

        // check edge lengths
        for (unsigned int i = 0; i < box.numEdges(); ++i)
            if ( abs(box.edge(i).length() - f) > eps)
                throw std::runtime_error("Edge " + std::to_string(i) + " has unexpected length");

        // check face areas
        for (unsigned int i = 0; i < box.numFaces(); ++i)
            if ( abs(box.face(i).area() - f*f) > eps*eps)
                throw std::runtime_error("Face " + std::to_string(i) + " has unexpected area");

        // check contains() query
        if (!box.contains(Point(0.5*f, 0.0, 0.0), eps))
            throw std::runtime_error("contains() query 1 failed");
        if (!box.contains(Point(0.5*f, 1e-6*f, 0.0), eps))
            throw std::runtime_error("contains() query 2 failed");
        if (box.contains(Point(0.5*f, -1e-6*f, 0.0), eps))
            throw std::runtime_error("contains() query 3 failed");
        if (!box.contains(Point(0.5*f, 0.5*f, f), eps))
            throw std::runtime_error("contains() query 4 failed");
        if (!box.contains(Point(0.5*f, 0.5*f, f-1e-6*f), eps))
            throw std::runtime_error("contains() query 5 failed");
        if (box.contains(Point(0.5*f, 0.5*f, f+1e-6*f), eps))
            throw std::runtime_error("contains() query 6 failed");

        // check correctness of center point
        const auto lowerCorner = box.corner(0);
        const auto upperCorner = box.corner(box.numCorners()-1);
        const auto diagonal = Segment(lowerCorner, upperCorner);
        if (!diagonal.center().isEqual(box.center(), eps))
            throw std::runtime_error("center point test failed");
    }

    std::cout << "All tests passed" << std::endl;
    return 0;
}
