#include <random>
#include <iostream>

#include <frackit/common/id.hh>
#include <frackit/io/gmshwriter.hh>

#include <frackit/sampling/makeuniformpointsampler.hh>
#include <frackit/sampling/quadrilateralsampler.hh>
#include <frackit/sampling/status.hh>

#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/box.hh>
#include <frackit/geometry/cylinder.hh>
#include <frackit/magnitude/containedmagnitude.hh>

#include <frackit/entitynetwork/constraints.hh>
#include <frackit/entitynetwork/networkbuilder.hh>

//! create a network of 3d quadrilaterals contained in a cylinder
int main()
{
    using namespace Frackit;

    // We use the same samplers for quadrilaterals as in example 1.
    // Therefore, we do not repeat the comments here. For more details
    // please see the main file of example1. We will only add comments
    // here in those places that deviate from the first example.
    static constexpr int worldDimension = 3;
    using ctype = double;
    using Quad = Quadrilateral<ctype, worldDimension>;

    // In this example we consider a cylindrical domain
    Cylinder<ctype> domain(/*radius*/0.5, /*height*/1.0);

    // However, we still sample the quadrilateral sample points
    // within the unit cube (a sampler on cylinders is also available)
    // Note that we shift the unit cube to have the origin in
    // the center of the bottom boundary
    Box<ctype> unitCube(-0.5, -0.5, 0.0, 0.5, 0.5, 1.0);

    using QuadSampler = QuadrilateralSampler<worldDimension>;
    using Distro = std::normal_distribution<ctype>;
    QuadSampler quadSampler1(makeUniformPointSampler(unitCube),
                             Distro(toRadians(45.0), toRadians(5.0)),
                             Distro(toRadians(90.0), toRadians(5.0)),
                             Distro(0.5, 0.1),
                             0.05);

    QuadSampler quadSampler2(makeUniformPointSampler(unitCube),
                             Distro(toRadians(0.0), toRadians(5.0)),
                             Distro(toRadians(0.0), toRadians(5.0)),
                             Distro(0.5, 0.1),
                             0.05);

    EntityNetworkConstraints<ctype> constraintsOnSelf;
    constraintsOnSelf.setMinDistance(0.05);
    constraintsOnSelf.setMinIntersectingAngle(toRadians(30.0));
    constraintsOnSelf.setMinIntersectionMagnitude(0.05);
    constraintsOnSelf.setMinIntersectionDistance(0.05);

    auto constraintsOnOther = constraintsOnSelf;
    constraintsOnOther.setMinIntersectingAngle(toRadians(40.0));

    // We can now also enforce constraints w.r.t. the domain boundary
    auto constraintsOnBoundary = constraintsOnSelf;

    std::vector<Quad> entitySet1;
    std::vector<Quad> entitySet2;
    const Id idSet1(1);
    const Id idSet2(2);

    SamplingStatus status;
    status.setTargetCount(idSet1, 12);
    status.setTargetCount(idSet2, 12);

    bool sampleIntoSet1 = true;

    std::cout << "\n --- Start entity sampling ---\n" << std::endl;
    while (!status.finished())
    {
        auto quad = sampleIntoSet1 ? quadSampler1() : quadSampler2();
        auto& entitySet = sampleIntoSet1 ? entitySet1 : entitySet2;
        const auto& otherEntitySet = sampleIntoSet1 ? entitySet2 : entitySet1;

        if (!constraintsOnSelf.evaluate(entitySet, quad))
        { status.increaseRejectedCounter(); continue; }

        if (!constraintsOnOther.evaluate(otherEntitySet, quad))
        { status.increaseRejectedCounter(); continue; }

        // we want to enforce constraints also w.r.t. to the cylinder boundary
        if (!constraintsOnBoundary.evaluate(domain.topFace(), quad))
        { status.increaseRejectedCounter(); continue; }
        if (!constraintsOnBoundary.evaluate(domain.bottomFace(), quad))
        { status.increaseRejectedCounter(); continue; }
        if (!constraintsOnBoundary.evaluate(OCCUtilities::getShape(domain.lateralFace()), quad))
        { status.increaseRejectedCounter(); continue; }

        // we also want to neglect quadrilaterals of which only
        // small fragments remain after confining them to the
        // domain. This computes the area of the quad that is
        // inside the domain
        const auto containedArea = computeContainedMagnitude(quad, domain);

        // reject if this is too small (< 0.01 m^2)
        if (containedArea < 0.01)
        { status.increaseRejectedCounter(); continue; }

        entitySet.push_back(quad);
        status.increaseCounter(sampleIntoSet1 ? idSet1 : idSet2);
        status.print();

        sampleIntoSet1 = !sampleIntoSet1;
    }
    std::cout << "\n --- Finished entity sampling ---\n" << std::endl;

    // We can now create a contained entity network from the two sets,
    // which has information on both the entities and the domain.
    ContainedEntityNetworkBuilder<ctype> builder;

    // define the domain (single sub-domain) and give it a unique id
    builder.addConfiningSubDomain(domain, Id(1));

    // define entities to be embedded in this domain
    builder.addSubDomainEntities(entitySet1, Id(1));
    builder.addSubDomainEntities(entitySet2, Id(1));

    const auto network = builder.build();

    // This can be written out in Gmsh (.geo) format to be
    // meshed by a three-dimensional mesh that is conforming
    // to the two-dimensional quadrilaterals
    std::cout << "\n --- Writing .geo file ---\n" << std::endl;
    GmshWriter writer(network);
    writer.write("network", // filename of the .geo files (will add extension .geo automatically)
                 0.1,       // element size to be used on the quadrilaterals
                 0.2);      // element size to be used in the domain away from the quadrilaterals
    std::cout << "\n --- Finished writing .geo file ---\n" << std::endl;

    return 0;
}
