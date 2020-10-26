#include <stdexcept>
#include <string>
#include <vector>
#include <random>
#include <cmath>

#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>

#include <frackit/common/math.hh>
#include <frackit/geometry/box.hh>
#include <frackit/occ/breputilities.hh>

#include <frackit/sampling/polygonsampler.hh>
#include <frackit/sampling/makeuniformpointsampler.hh>

//! test random sampling of polygons in 3d space
int main()
{
    using ctype = double;

    using namespace Frackit;
    using Box = Box<ctype>;
    using Polygon = Polygon<ctype, 3>;

    // Sample points within unit cube
    Box box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    PolygonSampler<3> polygonSampler(makeUniformPointSampler(box),
                                     std::normal_distribution<ctype>(toRadians(45.0), toRadians(0.05)),
                                     std::normal_distribution<ctype>(toRadians(45.0), toRadians(0.05)),
                                     std::uniform_real_distribution<ctype>(1.0, 2.0),
                                     std::uniform_real_distribution<ctype>(1.0, 2.0),
                                     std::uniform_int_distribution<int>(4, 10));

    // sample 10 polygons
    std::vector<Polygon> polygons;
    for (unsigned int i = 0; i < 10; ++i)
        polygons.emplace_back(polygonSampler());

    // create a single compound and write to .brep file
    // build a single compound shape
    BRep_Builder b;
    TopoDS_Compound c;
    b.MakeCompound(c);
    for (const auto& polygon : polygons)
        b.Add(c, OCCUtilities::getShape(polygon));
    b.Add(c, OCCUtilities::getShape(box));
    BRepTools::Write(c, "polygonsinbox.brep");

    std::cout << "All tests passed" << std::endl;
    return 0;
}
