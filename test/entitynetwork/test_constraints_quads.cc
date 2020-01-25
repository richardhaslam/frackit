#include <cmath>
#include <string>
#include <stdexcept>

#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/cylinder.hh>
#include <frackit/entitynetwork/constraints.hh>
#include <frackit/occ/breputilities.hh>

//! test the constraints for entity networks of quadrilaterals
int main()
{
    using ctype = double;

    using Quad = Frackit::Quadrilateral<ctype, 3>;
    using Point = typename Quad::Point;
    using Vector = typename Frackit::Vector<ctype, 3>;

    // Basis Vectors
    const Vector e1(1.0, 0.0, 0.0);
    const Vector e2(0.0, 1.0, 0.0);
    const Vector e3(0.0, 0.0, 1.0);

    // Define constraints
    Frackit::EntityNetworkConstraints<ctype> constraints;
    constraints.setMinDistance(0.1);
    constraints.setMinIntersectingAngle(M_PI/4.0);
    constraints.setMinIntersectionMagnitude(0.05);
    constraints.setMinIntersectionDistance(0.05);

    // disk to test the others against
    Quad mainQuad(Point(-1.0, -1.0, 0.0),
                  Point( 1.0, -1.0, 0.0),
                  Point(-1.0,  1.0, 0.0),
                  Point( 1.0,  1.0, 0.0));

    // violates distance constraint
    Quad quad1(Point(-1.0, -1.0, 0.1 - 1.0e-3),
               Point( 1.0, -1.0, 0.1 - 1.0e-3),
               Point(-1.0,  1.0, 0.1 - 1.0e-3),
               Point( 1.0,  1.0, 0.1 - 1.0e-3));
    if (constraints.evaluate(mainQuad, quad1))
        throw std::runtime_error("Did not detect distance violation");
    std::cout << "Test 1 passed" << std::endl;

    // just doesn't violates distance constraint
    Quad quad2(Point(-1.0, -1.0, 0.1 + 1.0e-3),
               Point( 1.0, -1.0, 0.1 + 1.0e-3),
               Point(-1.0,  1.0, 0.1 + 1.0e-3),
               Point( 1.0,  1.0, 0.1 + 1.0e-3));
    if (!constraints.evaluate(mainQuad, quad2))
        throw std::runtime_error("Detected false positive distance violation");
    std::cout << "Test 2 passed" << std::endl;

    // violates distance constraint (orthogonal)
    Quad quad3(Point(-1.0, 0.0, 0.1 - 1.0e-3),
               Point( 1.0, 0.0, 0.1 - 1.0e-3),
               Point(-1.0, 0.0, 1.1 - 1.0e-3),
               Point( 1.0, 0.0, 1.1 - 1.0e-3));
    if (constraints.evaluate(mainQuad, quad3))
        throw std::runtime_error("Did not detect distance violation");
    std::cout << "Test 3 passed" << std::endl;

    // just doesn't violates distance constraint (orthogonal)
    Quad quad4(Point(-1.0, 0.0, 0.1 + 1.0e-3),
               Point( 1.0, 0.0, 0.1 + 1.0e-3),
               Point(-1.0, 0.0, 1.1 + 1.0e-3),
               Point( 1.0, 0.0, 1.1 + 1.0e-3));
    if (!constraints.evaluate(mainQuad, quad4))
        throw std::runtime_error("Detected false positive distance violation");
    std::cout << "Test 4 passed" << std::endl;

    // too small intersection
    Quad quad5(Point(-0.05+1e-3, 0.0,  1.0),
               Point( 0.0     , 0.0, -1.0),
               Point( 0.0     , 0.0,  1.5),
               Point( 0.05-1e-3, 0.0,  1.0));
    if (constraints.evaluate(mainQuad, quad5))
        throw std::runtime_error("Did not detect intersection magnitude violation");
    std::cout << "Test 5 passed" << std::endl;

    // intersection magnitude ok
    Quad quad6(Point(-0.05-1e-3, 0.0,  1.0),
               Point( 0.0     , 0.0, -1.0),
               Point( 0.0     , 0.0,  1.5),
               Point( 0.05+1e-3, 0.0,  1.0));
    if (!constraints.evaluate(mainQuad, quad6))
        throw std::runtime_error("False positive intersection magnitude violation");
    std::cout << "Test 6 passed" << std::endl;

    // intersection magnitude ok, but distance to boundary is too small
    Quad quad7(Point(-2.0 + 0.0999, 0.0,  1.0),
               Point( 0.0,          0.0, -1.0),
               Point( 0.0,          0.0,  1.5),
               Point( 2.0 - 0.0999, 0.0,  1.0));
    if (constraints.evaluate(mainQuad, quad7))
        throw std::runtime_error("Did not detect intersection distance violation");
    std::cout << "Test 7 passed" << std::endl;

    // intersection magnitude ok, distance to boundary is just ok
    Quad quad8(Point(-2.0 + 0.1001, 0.0,  1.0),
               Point( 0.0,          0.0, -1.0),
               Point( 0.0,          0.0,  1.5),
               Point( 2.0 - 0.1001, 0.0,  1.0));
    if (!constraints.evaluate(mainQuad, quad8))
        throw std::runtime_error("False positive intersection distance violation");
    std::cout << "Test 8 passed" << std::endl;

    // test the two above also with a shape representation of the main quad
    {
        if (constraints.evaluate(mainQuad, Frackit::OCCUtilities::getShape(quad7)))
            throw std::runtime_error("Did not detect intersection distance violation");
        std::cout << "Test 7 (with shape) passed" << std::endl;

        if (!constraints.evaluate(mainQuad, Frackit::OCCUtilities::getShape(quad8)))
            throw std::runtime_error("False positive intersection distance violation");
        std::cout << "Test 8 (with shape) passed" << std::endl;
    }

    // intersection angle just ok
    Quad quad9(Point(-0.5, -0.5, -0.5 - 1e-3),
               Point( 0.5, -0.5, -0.5 - 1e-3),
               Point(-0.5,  0.5,  0.5 + 1e-3),
               Point( 0.5,  0.5,  0.5 + 1e-3));
    if (!constraints.evaluate(mainQuad, quad9))
        throw std::runtime_error("False positive intersection angle violation");
    std::cout << "Test 9 passed" << std::endl;

    // intersection angle too small
    Quad quad10(Point(-0.5, -0.5, -0.5 + 1e-3),
                Point( 0.5, -0.5, -0.5 + 1e-3),
                Point(-0.5,  0.5,  0.5 - 1e-3),
                Point( 0.5,  0.5,  0.5 - 1e-3));
    if (constraints.evaluate(mainQuad, quad10))
        throw std::runtime_error("Did not detect intersection angle violation");
    std::cout << "Test 10 passed" << std::endl;

    // Test constraints w.r.t. cylinder
    Frackit::Cylinder<ctype> cylinder(0.5, 1.0);
    Quad quad11(Point(-0.6, -0.0, 1.0 - 0.05111),
                Point( 0.0, -0.6, 1.0 - 0.05111),
                Point( 0.6,  0.0, 1.0 - 0.05111),
                Point( 0.0,  0.6, 1.0 - 0.05111));
    if (!constraints.evaluate(cylinder.lateralFace(), quad11))
        throw std::runtime_error("False positive intersection distance violation");
    std::cout << "Test 11 passed" << std::endl;

    // violates intersection distance constraint
    Quad quad12(Point(-0.6, -0.0, 1.0 - 0.04999),
                Point( 0.0, -0.6, 1.0 - 0.04999),
                Point( 0.6,  0.0, 1.0 - 0.04999),
                Point( 0.0,  0.6, 1.0 - 0.04999));
    if (constraints.evaluate(cylinder.lateralFace(), quad12))
        throw std::runtime_error("Did not detect intersection distance violation");
    std::cout << "Test 12 passed" << std::endl;

    // just doesn't violate distance constraint
    Quad quad13(Point(-0.40001, -0.0,     01.5),
                Point( 0.0,     -0.40001, 01.5),
                Point( 0.40001,  0.0,     01.5),
                Point( 0.0,      0.40001, 01.5));
    if (!constraints.evaluate(cylinder.lateralFace(), quad13))
        throw std::runtime_error("False positive distance violation");
    std::cout << "Test 13 passed" << std::endl;

    // just violates distance constraint
    Quad quad14(Point(-0.39999, -0.0,     0.5),
                Point( 0.0,     -0.39999, 0.5),
                Point( 0.39999,  0.0,     0.5),
                Point( 0.0,      0.39999, 0.5));
    if (constraints.evaluate(cylinder.lateralFace(), quad14))
        throw std::runtime_error("Did not detect distance violation");
    std::cout << "Test 14 passed" << std::endl;

    std::cout << "All tests passed" << std::endl;
    return 0;
}
