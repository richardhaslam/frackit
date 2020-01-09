#include <stdexcept>
#include <string>
#include <array>

#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>

#include <frackit/sampling/pointsampling.hh>
#include <frackit/geometry/cylinder.hh>
#include <frackit/occ/breputilities.hh>

//! test random sampling of points on geometries
int main()
{
    using ctype = double;

    using namespace Frackit;
    using Cylinder = Cylinder<ctype>;
    using CylPointSampler = GeometryPointSampler<Cylinder>;

    Cylinder cylinder(0.5, 1.0);
    CylPointSampler cylPointSampler(cylinder);

    // sample 500 points
    std::array<typename CylPointSampler::Point, 500> points;
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

    std::cout << "All tests passed" << std::endl;
    return 0;
}
