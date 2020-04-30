#include <stdexcept>
#include <string>
#include <array>

#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>

#include <frackit/sampling/makeuniformpointsampler.hh>
#include <frackit/geometry/point.hh>
#include <frackit/geometry/cylinder.hh>
#include <frackit/geometry/hollowcylinder.hh>
#include <frackit/occ/breputilities.hh>

//! test random sampling of points on geometries
int main()
{
    using ctype = double;

    using namespace Frackit;
    using Cylinder = Cylinder<ctype>;

    Cylinder cylinder(0.5, 1.0);
    auto cylPointSampler = makeUniformPointSampler(cylinder);

    // sample 500 points
    std::array<Frackit::Point<ctype, 3>, 500> points;
    for (unsigned int i = 0; i < 500; ++i)
        points[i] = cylPointSampler();

    // create a single compound and write to .brep file
    // build a single compound shape
    BRep_Builder b;
    TopoDS_Compound c;
    b.MakeCompound(c);
    for (const auto& p : points)
        b.Add(c, OCCUtilities::getShape(p));
    b.Add(c, OCCUtilities::getShape(cylinder));

    BRepTools::Write(c, "cylinderpoints.brep");

    // test sampling on a hollow cylinder
    HollowCylinder<ctype> hollowCylinder(0.4, 0.5, 1.0);
    cylPointSampler = makeUniformPointSampler(hollowCylinder);
    for (unsigned int i = 0; i < 500; ++i)
        points[i] = cylPointSampler();

    // create a single compound and write to .brep file
    // build a single compound shape
    TopoDS_Compound c2;
    b.MakeCompound(c2);
    for (const auto& p : points)
        b.Add(c2, OCCUtilities::getShape(p));
    b.Add(c2, OCCUtilities::getShape(cylinder));
    BRepTools::Write(c2, "cylinderpoints_hollow.brep");

    std::cout << "All tests passed" << std::endl;
    return 0;
}
