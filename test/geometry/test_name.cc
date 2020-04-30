#include <stdexcept>
#include <string>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/box.hh>
#include <frackit/geometry/vector.hh>
#include <frackit/geometry/line.hh>
#include <frackit/geometry/triangle.hh>

#include <frackit/geometry/cylinder.hh>
#include <frackit/geometry/ellipsearc.hh>
#include <frackit/geometry/circle.hh>

//! test virtual name interface
int main()
{
    using namespace Frackit;
    using ctype = double;

    // Except for triangles, we can test
    // all linear geometry types by extracting
    // from a box. This saves us from having to
    // instantiate all of them manually here.
    using Box = Box<ctype>;
    Box box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);

    // Test virtual name interface by casting all
    // geometry types into pointer on abstract base class
    Geometry* geometry = &box;
    if (geometry->name() != "Box") throw std::runtime_error("Box name wrong");

    // try casting geometry pointer back into box, get a point and set pointer
    Box* castBox = dynamic_cast<Box*>(geometry);
    auto point = castBox->corner(1);
    geometry = &point;
    if (geometry->name() != "Point_3d") throw std::runtime_error("Point name wrong");

    // test segments
    auto segment = box.edge(1);
    geometry = &segment;
    if (geometry->name() != "Segment_3d") throw std::runtime_error("Segment name wrong");

    // test line
    auto line = segment.supportingLine();
    geometry = &line;
    if (geometry->name() != "Line_3d") throw std::runtime_error("Line name wrong");

    // test direction
    Line<ctype, 3>* lineCast = dynamic_cast<Line<ctype, 3>*>(geometry);
    auto direction = lineCast->direction();
    geometry = &direction;
    if (geometry->name() != "Direction_3d") throw std::runtime_error("Direction name wrong");

    // test vector
    auto vector = Vector<ctype, 3>(direction);
    geometry = &vector;
    if (geometry->name() != "Vector_3d") throw std::runtime_error("Vector name wrong");

    // test quadrilateral
    auto face = box.face(1);
    geometry = &face;
    if (geometry->name() != "Quadrilateral_3d") throw std::runtime_error("Quadrilateral name wrong");

    // test plane
    auto plane = face.supportingPlane();
    geometry = &plane;
    if (geometry->name() != "Plane_3d") throw std::runtime_error("Plane name wrong");

    // test triangle
    Triangle<ctype, 3> triangle(face.corner(0), face.corner(1), face.corner(2));
    geometry = &triangle;
    if (geometry->name() != "Triangle_3d") throw std::runtime_error("Triangle name wrong");

    // Test the non-linear geometries by extracting them from a cylinder
    Cylinder<ctype> cylinder(0.5, 1.0);
    geometry = &cylinder;
    if (geometry->name() != "Cylinder") throw std::runtime_error("Cylinder name wrong");

    // test disk
    auto disk = cylinder.bottomFace();
    geometry = &disk;
    if (geometry->name() != "Disk") throw std::runtime_error("Disk name wrong");

    // test ellipse
    auto ellipse = disk.boundingEllipse();
    geometry = &ellipse;
    if (geometry->name() != "Ellipse_3d") throw std::runtime_error("Ellipse name wrong");

    // test ellipse arc
    auto p1 = ellipse.getPoint(0.0);
    auto p2 = ellipse.getPoint(0.25);
    EllipseArc<ctype, 3> arc(ellipse, p1, p2);
    geometry = &arc;
    if (geometry->name() != "EllipseArc_3d") throw std::runtime_error("EllipseArc name wrong");

    // test cylinder surface
    auto cylMantle = cylinder.lateralFace();
    geometry = &cylMantle;
    if (geometry->name() != "CylinderMantle") throw std::runtime_error("CylinderMantle name wrong");

    // test circle
    Circle<ctype, 3> circle(p1, plane.normal(), 1.0);
    geometry = &circle;
    if (geometry->name() != "Circle_3d") throw std::runtime_error("Circle name wrong");

    std::cout << "All tests passed" << std::endl;
    return 0;
}
