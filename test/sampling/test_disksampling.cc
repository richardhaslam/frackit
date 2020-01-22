#include <stdexcept>
#include <string>
#include <vector>
#include <random>

#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>

#include <frackit/common/math.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/occ/breputilities.hh>

#include <frackit/sampling/disksampler.hh>
#include <frackit/sampling/makeuniformpointsampler.hh>

//! test random sampling of points on geometries
int main()
{
    using ctype = double;

    using namespace Frackit;
    using Cylinder = Cylinder<ctype>;
    using Disk = Disk<ctype>;

    Cylinder cylinder(0.5, 1.0);
    auto cylPointSampler = makeUniformPointSampler(cylinder);
    DiskSampler diskSampler(std::normal_distribution<ctype>(0.5, 0.05),
                            std::normal_distribution<ctype>(0.25, 0.025),
                            std::normal_distribution<ctype>(toRadians(45.0),toRadians(10.0)),
                            std::normal_distribution<ctype>(toRadians(45.0),toRadians(10.0)),
                            std::normal_distribution<ctype>(toRadians(45.0),toRadians(10.0)));

    // sample 5 disks
    std::vector<Disk> disks;
    for (unsigned int i = 0; i < 5; ++i)
        disks.emplace_back(diskSampler( cylPointSampler() ));

    // create a single compound and write to .brep file
    // build a single compound shape
    BRep_Builder b;
    TopoDS_Compound c;
    b.MakeCompound(c);
    for (const auto& disk : disks)
        b.Add(c, OCCUtilities::getShape(disk));
    b.Add(c, OCCUtilities::getShape(cylinder));

    BRepTools::Write(c, "disksincylinder.brep");
    std::cout << "All tests passed" << std::endl;
    return 0;
}
