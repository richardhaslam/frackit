#include <random>
#include <iostream>

#include <frackit/common/id.hh>
#include <frackit/io/gmshwriter.hh>

#include <frackit/sampling/makeuniformpointsampler.hh>
#include <frackit/sampling/quadrilateralsampler.hh>
#include <frackit/sampling/status.hh>

#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/box.hh>

#include <frackit/entitynetwork/constraints.hh>
#include <frackit/entitynetwork/networkbuilder.hh>

//! create a network of 3d quadrilaterals
int main()
{
    using namespace Frackit;

    // We consider 3d space here
    static constexpr int worldDimension = 3;

    // Define the type used for coordinates
    using ctype = double;

    // Internal geometry type for 3d quadrilaterals
    using Quad = Quadrilateral<ctype, worldDimension>;

    // Define a domain (here: unit cube) in which the quadrilaterals should be created.
    // Boxes are created by providing xmin, ymin, zmin and xmax, ymax and zmax in constructor.
    Box<ctype> domain(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);

    // We now create a sampler instance that uniformly samples points within this box.
    // These points will be used as the quadrilateral centers.
    auto pointSampler = makeUniformPointSampler(domain);

    // Sampler class for quadrilaterals. Per default, this uses
    // uniform distributions for all parameters defining the quadrilaterals.
    // Quadrilateral samplers require distributions for strike angle, dip angle,
    // edge length (see class description for more details). Moreover, we define
    // a minimum edge length.
    using QuadSampler = QuadrilateralSampler<worldDimension>;
    using Distro = std::normal_distribution<ctype>;
    QuadSampler quadSampler1(pointSampler,
                             Distro(toRadians(45.0), toRadians(5.0)),  // strike angle: mean value & standard deviation
                             Distro(toRadians(90.0), toRadians(5.0)),  // dip angle: mean value & standard deviation
                             Distro(0.5, 0.1),                         // edge length: mean value & standard deviation
                             0.05);                                    // threshold for minimum edge length

    // We use this sampler to create quadrilaterals based on the distributions.
    // However, we want to create another set of quadrilaterals, whose entities
    // are approximately orthogonal to those defined with the above sampler.
    QuadSampler quadSampler2(makeUniformPointSampler(domain),         // use a new point sampler instance!
                             Distro(toRadians(0.0), toRadians(5.0)),  // strike angle: mean value & standard deviation
                             Distro(toRadians(0.0), toRadians(5.0)),  // dip angle: mean value & standard deviation
                             Distro(0.5, 0.1),                        // edge length: mean value & standard deviation
                             0.05);                                   // threshold for minimum edge length

    // We want to enforce some constraints on the set of quadrilaterals.
    // In particular, for entities of the same set we want a minimum spacing
    // distance of 5cm, and the quadrilaterals must not intersect in angles
    // less than 30°. Moreover, if they intersect, we don't want intersection
    // edges whose length is smaller than 5cm, and, the intersection should not
    // be too close to the boundary of one of two intersecting quadrilaterals. Here: 5cm.
    EntityNetworkConstraints constraintsOnSelf;
    constraintsOnSelf.setMinDistance(0.05);
    constraintsOnSelf.setMinIntersectingAngle(toRadians(30.0));
    constraintsOnSelf.setMinIntersectionMagnitude(0.05);
    constraintsOnSelf.setMinIntersectionDistance(0.05);

    // with respect to entities of the other set, we want to have larger intersection angles
    auto constraintsOnOther = constraintsOnSelf;
    constraintsOnOther.setMinIntersectingAngle(toRadians(40.0));

    // container to store created entities
    std::vector<Quad> entitySet1;
    std::vector<Quad> entitySet2;

    // we give unique identifiers to both entity sets
    const Id idSet1(1);
    const Id idSet2(2);

    // use the status class to define when to stop sampling
    SamplingStatus status;
    status.setTargetCount(idSet1, 15); // we want 15 entities in set 1
    status.setTargetCount(idSet2, 15); // we want 15 entities in set 2

    // start sampling into set 1 and keep alternating
    bool sampleIntoSet1 = true;

    std::cout << "\n --- Start entity sampling ---\n" << std::endl;
    while (!status.finished())
    {
        // sample a quadrilateral, alternating between sampler 1 and sampler 2
        auto quad = sampleIntoSet1 ? quadSampler1() : quadSampler2();
        auto& entitySet = sampleIntoSet1 ? entitySet1 : entitySet2;
        const auto& otherEntitySet = sampleIntoSet1 ? entitySet2 : entitySet1;

        // sample again if constraints w.r.t. other
        // entities of this set are not fulfilled
        if (!constraintsOnSelf.evaluate(entitySet, quad))
        { status.increaseRejectedCounter(); continue; }

        // sample again if constraints w.r.t. other
        // entities of the other set are not fulfilled
        if (!constraintsOnOther.evaluate(otherEntitySet, quad))
        { status.increaseRejectedCounter(); continue; }

        // the quadrilateral is admissible
        entitySet.push_back(quad);

        // tell the status we have a new entity in this set
        const auto& setId = sampleIntoSet1 ? idSet1 : idSet2;
        status.increaseCounter(setId);
        status.print();

        // sample into the other set the next time
        sampleIntoSet1 = !sampleIntoSet1;
    }
    std::cout << "\n --- Finished entity sampling ---\n" << std::endl;

    // We can now create an entity network from the two sets
    EntityNetworkBuilder builder;
    builder.addEntities(entitySet1);
    builder.addEntities(entitySet2);

    const auto network = builder.build();

    // This can be written out in Gmsh (.geo) format to be
    // meshed by a two-dimensional surface mesh
    std::cout << "\n --- Writing .geo file ---\n" << std::endl;
    GmshWriter writer(network);
    writer.write("network", // filename of the .geo files (will add extension .geo automatically)
                 0.1);      // element size to be used
    std::cout << "\n --- Finished writing .geo file ---\n" << std::endl;

    return 0;
}
