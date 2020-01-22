#include <iostream>
#include <random>
#include <fstream>
#include <unordered_map>

#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_DataMapOfShapeInteger.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS.hxx>

#include <frackit/geometry/cylinder.hh>
#include <frackit/geometry/disk.hh>

#include <frackit/magnitude/magnitude.hh>
#include <frackit/magnitude/containedmagnitude.hh>

#include <frackit/common/id.hh>
#include <frackit/occ/breputilities.hh>
#include <frackit/precision/precision.hh>

#include <frackit/sampling/makeuniformpointsampler.hh>
#include <frackit/sampling/disksampler.hh>

#include <frackit/entitynetwork/containedentitynetwork.hh>
#include <frackit/entitynetwork/networkbuilder.hh>
#include <frackit/io/gmshwriter.hh>

//! create 3 disk-shaped fractures and write geo (gmsh file format) file
int main(int argc, char** argv)
{
    using namespace Frackit;
    using ctype = double;
    using Domain = Cylinder<ctype>;
    using Disk = Disk<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;

    const Direction e1(Vector(1.0, 0.0, 0.0));
    const Direction e2(Vector(0.0, 1.0, 0.0));
    const Direction e3(Vector(0.0, 0.0, 1.0));

    // Domain in which the disks are embedded
    Domain domain(0.5, 1.0);
    Domain domain2(Disk(Point(0.0, 0.0, 1.0), e1, e2, 1.0, 1.0), 1.0);

    ContainedEntityNetworkBuilder builder;
    builder.addConfiningSubDomain(domain, Id(1));
    builder.addSubDomainEntity(Disk(Point(0.0, 0.0, 0.1), e1, e2, 1.0, 1.0),  Id(1));
    builder.addSubDomainEntity(Disk(Point(0.0, 0.0, 0.5), e1, e3, 2.0, 2.0),  Id(1));
    builder.addSubDomainEntity(Disk(Point(0.0, 0.0, 0.75), e1, e2, 1.0, 1.0), Id(1));

    builder.addConfiningSubDomain(domain2, Id(2));
    builder.addSubDomainEntity(Disk(Point(0.0, 0.0, 1.1), e1, e2, 1.0, 1.0),  Id(2));
    builder.addSubDomainEntity(Disk(Point(0.0, 0.0, 1.5), e1, e3, 2.0, 2.0),  Id(2));
    builder.addSubDomainEntity(Disk(Point(0.0, 0.0, 1.75), e1, e2, 1.0, 1.0), Id(2));

    GmshWriter writer(builder.build());
    writer.write("final", 1.0, 1.0);

    std::cout << "All tests passed" << std::endl;

    return 0;
}
