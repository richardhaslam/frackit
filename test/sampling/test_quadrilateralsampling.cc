#include <stdexcept>
#include <string>
#include <vector>
#include <random>

#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>

#include <frackit/common/math.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/box.hh>
#include <frackit/occ/breputilities.hh>

#include <frackit/sampling/quadrilateralsampler.hh>
#include <frackit/sampling/makeuniformpointsampler.hh>

//! test random sampling of quadrilaterals in 3d space
int main()
{
    using ctype = double;

    using namespace Frackit;
    using Box = Box<ctype>;
    using Quad = Quadrilateral<ctype, 3>;

    // Sample points within unit cube
    Box box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);

    auto pointSampler = makeUniformPointSampler(box);
    QuadrilateralSampler<3> quadSampler(makeUniformPointSampler(box),
                                        std::normal_distribution<ctype>(toRadians(45.0), toRadians(5.0)),
                                        std::normal_distribution<ctype>(toRadians(20.0), toRadians(5.0)),
                                        std::uniform_real_distribution<ctype>(0.4, 0.6),
                                        std::uniform_real_distribution<ctype>(0.4, 0.6));

    // sample 10 quadrilaterals
    std::vector<Quad> quads;
    for (unsigned int i = 0; i < 5; ++i)
        quads.emplace_back(quadSampler());

    // create a single compound and write to .brep file
    // build a single compound shape
    BRep_Builder b;
    TopoDS_Compound c;
    b.MakeCompound(c);
    for (const auto& quad : quads)
        b.Add(c, OCCUtilities::getShape(quad));
    b.Add(c, OCCUtilities::getShape(box));

    BRepTools::Write(c, "quadsinbox.brep");
    std::cout << "All tests passed" << std::endl;
    return 0;
}
