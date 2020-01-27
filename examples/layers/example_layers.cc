#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <random>

// basic geometry types
#include <frackit/geometry/box.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/quadrilateral.hh>

// headers containing functions to compute lengths/areas/volumes
#include <frackit/magnitude/magnitude.hh>
#include <frackit/magnitude/containedmagnitude.hh>

// sampler for points and disks
#include <frackit/sampling/makeuniformpointsampler.hh>
#include <frackit/sampling/disksampler.hh>
#include <frackit/sampling/quadrilateralsampler.hh>
#include <frackit/sampling/multigeometrysampler.hh>
#include <frackit/sampling/sequentialsamplingstrategy.hh>
#include <frackit/sampling/status.hh>

// constraints to be enforced on the network (distance, angles, etc.)
#include <frackit/entitynetwork/constraints.hh>
#include <frackit/entitynetwork/constraintsmatrix.hh>
#include <frackit/entitynetwork/multigeometryentityset.hh>

// builder class for creating networks of entities confined in (sub-)domains
#include <frackit/entitynetwork/networkbuilder.hh>

// writes an entity network to a meshable Gmsh .geo file format
#include <frackit/io/gmshwriter.hh>

int main(int argc, char** argv)
{
    //! print welcome message
    std::cout << "\n\n"
              << "##############################################################################\n"
              << "## Example: Creation of entities (disks/quadrilaterals) in a layered medium ##\n"
              << "##############################################################################\n"
              << "\n\n";

    // Define some types used here
    using namespace Frackit;
    using ctype = double;
    using Disk = Disk<ctype>;
    using Quad = Quadrilateral<ctype, 3>;

    /////////////////////////////////////////////////////
    // 1. Read in the domain geometry from .brep file. //
    //    The file name is defined in CMakeLists.txt   //
    /////////////////////////////////////////////////////
    const auto domainShape = OCCUtilities::readShape(BREPFILE);

    // obtain the three solids contained in the file
    const auto solids = OCCUtilities::getSolids(domainShape);
    if (solids.size() != 3)
        throw std::runtime_error("Expected the .brep file to contain 3 solids");

    // The sub-domain we want to create a network in is the center one.
    // Compute its volume and get the boundary faces for constraint evaluation.
    const auto& networkDomain = solids[1];
    const auto domainVolume = computeMagnitude(networkDomain);
    const auto shells = OCCUtilities::getShells(networkDomain);
    if (shells.size() != 1) throw std::runtime_error("Expected a single shell bounding the domain");
    const auto domainBoundaryFaces = OCCUtilities::getFaces(shells[0]);

    //////////////////////////////////////////////////////////////////////////////
    // 2. Make sampler classes to randomly sample points (entity centers)       //
    //    and disk/quadrilaterals as entities (using the sampled center points) //
    //////////////////////////////////////////////////////////////////////////////

    // Bounding box of the domain in which we want to place the entities
    const auto domainBBox = OCCUtilities::getBoundingBox(networkDomain);

    // sampler for disks (orientation 1)
    DiskSampler diskSampler(makeUniformPointSampler(domainBBox),                               // sampler for disk center points
                            std::normal_distribution<ctype>(20.0, 6.5),                        // major axis length: mean value & standard deviation
                            std::normal_distribution<ctype>(15.0, 4.5),                        // minor axis length: mean value & standard deviation
                            std::normal_distribution<ctype>(toRadians(0.0), toRadians(7.5)),   // rotation around x-axis: mean value & standard deviation
                            std::normal_distribution<ctype>(toRadians(0.0), toRadians(7.5)),   // rotation around y-axis: mean value & standard deviation
                            std::normal_distribution<ctype>(toRadians(0.0), toRadians(7.5)));  // rotation around z-axis: mean value & standard deviation

    // sampler for quadrilaterals (orientation 2)
    QuadrilateralSampler<3> quadSampler(makeUniformPointSampler(domainBBox),                               // sampler for quadrilateral center points
                                        std::normal_distribution<ctype>(toRadians(0.0), toRadians(5.0)),   // strike angle: mean value & standard deviation
                                        std::normal_distribution<ctype>(toRadians(45.0), toRadians(5.0)),  // dip angle: mean value & standard deviation
                                        std::normal_distribution<ctype>(40.0, 5.0),                        // edge length: mean value & standard deviation
                                        5.0);                                                              // threshold for minimum edge length

    // Define ids for the two entity sets
    const Id diskSetId(1); // we give the set of orientation one, consisting of disks, the id 1
    const Id quadSetId(2); // we give the set of orientation two, consisting of quadrilaterals, the id 2

    // Sampler that samples both geometry types (disks & quads)
    // In this one can define an arbitrary number of samplers, each
    // of which is associated with an entity set with a unique id
    MultiGeometrySampler<Disk, Quad> multiSampler;
    multiSampler.addGeometrySampler(diskSampler, diskSetId);
    multiSampler.addGeometrySampler(quadSampler, quadSetId);

    // Define constraints between entities of the same set
    using Constraints = EntityNetworkConstraints<ctype>;
    Constraints constraintsOnSelf;
    constraintsOnSelf.setMinDistance(2.5);
    constraintsOnSelf.setMinIntersectingAngle(toRadians(25.0));
    constraintsOnSelf.setMinIntersectionMagnitude(2.5);
    constraintsOnSelf.setMinIntersectionDistance(2.0);

    // Define constraints between entities of different sets
    Constraints constraintsOnOther;
    constraintsOnOther.setMinDistance(2.5);
    constraintsOnOther.setMinIntersectingAngle(toRadians(40.0));
    constraintsOnOther.setMinIntersectionMagnitude(2.5);
    constraintsOnOther.setMinIntersectionDistance(2.0);

    // We can use a constraint matrix to facilitate constraint evaluation
    EntityNetworkConstraintsMatrix<Constraints> constraintsMatrix;
    constraintsMatrix.addConstraints(constraintsOnSelf,               // constraint instance
                                     {IdPair(diskSetId, diskSetId),   // sets between which to use these constraints
                                      IdPair(quadSetId, quadSetId)}); // sets between which to use these constraints

    constraintsMatrix.addConstraints(constraintsOnOther,              // constraint instance
                                     {IdPair(diskSetId, quadSetId),   // sets between which to use these constraints
                                      IdPair(quadSetId, diskSetId)}); // sets between which to use these constraints

    // Moreover, we define constraints w.r.t. the domain boundary
    EntityNetworkConstraints constraintsOnDomain;
    constraintsOnDomain.setMinDistance(1.5);
    constraintsOnDomain.setMinIntersectingAngle(toRadians(5.0));
    constraintsOnDomain.setMinIntersectionMagnitude(20.0);
    constraintsOnDomain.setMinIntersectionDistance(2.0);

    // Helper class to store an arbitrary number
    // of entity sets of different geometry types.
    MultiGeometryEntitySet<Disk, Quad> entitySets;

    // Helper class for terminal output of the creation
    // progress and definition of stop criterion etc
    SamplingStatus status;
    status.setTargetCount(diskSetId, 15); // we want 15 entities of orientation 1
    status.setTargetCount(quadSetId, 15); // we want 15 entities of orientation 2

    ctype containedNetworkArea = 0.0;
    while (!status.finished())
    {
        // randomly sample a geometry and obtain the id
        // of the entity set for which it was sampled.
        Id id;
        auto geom = multiSampler(id);

        // If the set this geometry belongs to is finished, skip the rest
        if (status.finished(id))
        { status.increaseRejectedCounter(); continue; }

        // Moreover, we want to avoid small fragments (< 200 m²)
        const auto containedArea = computeContainedMagnitude(geom, networkDomain);
        if (containedArea < 200.0)
        { status.increaseRejectedCounter(); continue; }

        // enforce constraints w.r.t. to the other entities
        if (!constraintsMatrix.evaluate(entitySets, geom, id))
        { status.increaseRejectedCounter(); continue; }

        // enforce constraints w.r.t. the domain boundaries
        if (!constraintsOnDomain.evaluate(domainBoundaryFaces, geom))
        { status.increaseRejectedCounter(); continue; }

        // the disk is admissible
        entitySets.addEntity(geom, id);
        status.increaseCounter(id);
        status.print();

        // keep track of entity area
        containedNetworkArea += containedArea;
    }

    // compute entity density
    const auto density = containedNetworkArea/domainVolume;
    std::cout << "\nEntity density of the contained network: " << density << " m²/m³" << std::endl;

    //////////////////////////////////////////////////////////////////////////
    // 3. The entities of the network have been created. We can now         //
    //    construct different types of networks (contained, confined, etc.) //
    //////////////////////////////////////////////////////////////////////////

    ContainedEntityNetworkBuilder containedConfinedBuilder;
    entitySets.exportEntitySets({diskSetId, quadSetId}, containedConfinedBuilder, Id(2));
    containedConfinedBuilder.addConfiningSubDomain(solids[0],     Id(1));
    containedConfinedBuilder.addConfiningSubDomain(networkDomain, Id(2));
    containedConfinedBuilder.addConfiningSubDomain(solids[2],     Id(3));

    // build network
    const auto containedConfinedNetwork = containedConfinedBuilder.build();

    // std::cout << "Constructing contained, confined network" << std::endl;
    // // 3.1 Construct a network such that the entities are confined in the
    // //     sub-domain and which contains information about the containing domain.
    // ContainedEntityNetworkBuilder containedConfinedBuilder;
    //
    // // define sub-domains
    // containedConfinedBuilder.addConfiningSubDomain(solids[0],     Id(1));
    // containedConfinedBuilder.addConfiningSubDomain(networkDomain, Id(2));
    // containedConfinedBuilder.addConfiningSubDomain(solids[2],     Id(3));
    //
    // // define entity network for sub-domain 2
    // containedConfinedBuilder.addSubDomainEntities(diskSet, Id(2));
    // containedConfinedBuilder.addSubDomainEntities(quadSet, Id(2));
    //
    // // build network
    // const auto containedConfinedNetwork = containedConfinedBuilder.build();
    //
    // std::cout << "Constructing contained, unconfined network" << std::endl;
    // // 3.2 Construct a network such that the entities are NOT confined in the
    // //     sub-domain and which contains information about the containing domain.
    // ContainedEntityNetworkBuilder containedUnconfinedBuilder;
    //
    // // define sub-domains
    // containedUnconfinedBuilder.addSubDomain(solids[0],     Id(1));
    // containedUnconfinedBuilder.addSubDomain(solids[2],     Id(3));
    // containedUnconfinedBuilder.addSubDomain(networkDomain, Id(2));
    //
    // // define entity network for sub-domain 2
    // containedUnconfinedBuilder.addSubDomainEntities(diskSet, Id(2));
    // containedUnconfinedBuilder.addSubDomainEntities(quadSet, Id(2));
    //
    // // build network
    // const auto containedUnconfinedNetwork = containedUnconfinedBuilder.build();
    //
    // std::cout << "Constructing uncontained, confined network" << std::endl;
    // // 3.3 Construct a network such that the entities are confined in the
    // //     sub-domain, but, which does not contain information about the domains.
    // //     This is computationally more efficient if only the 2d network is of interest.
    // //     Moreover, the files written to disk are smaller as the domain is not contained.
    // EntityNetworkBuilder uncontainedConfinedBuilder;
    //
    // // define sub-domains
    // uncontainedConfinedBuilder.addConfiningSubDomain(networkDomain, Id(1));
    //
    // // define entity network for sub-domain 2
    // uncontainedConfinedBuilder.addSubDomainEntities(diskSet, Id(1));
    // uncontainedConfinedBuilder.addSubDomainEntities(quadSet, Id(1));
    //
    // // build network
    // const auto uncontainedConfinedNetwork = uncontainedConfinedBuilder.build();
    //
    // std::cout << "Constructing uncontained, unconfined network" << std::endl;
    // // 3.4 Construct the network only, independent of any domains
    // EntityNetworkBuilder uncontainedUnconfinedBuilder;
    // uncontainedUnconfinedBuilder.addEntities(diskSet);
    // uncontainedUnconfinedBuilder.addEntities(quadSet);
    //
    // // build network
    // const auto uncontainedUnconfinedNetwork = uncontainedUnconfinedBuilder.build();
    //
    //
    // ///////////////////////////////////////////////////////
    // // 4. Write the networks to disk in Gmsh .geo format //
    // ///////////////////////////////////////////////////////
    std::cout << "Write networks to disk" << std::endl;
    GmshWriter containedConfinedWriter(containedConfinedNetwork);
    containedConfinedWriter.write("contained_confined", // body of the filename to be used (will add .geo)
                                  2.5,                  // mesh size to be used on entities
                                  5.0);                 // mesh size to be used on domain boundaries
    //
    // GmshWriter containedUnconfinedWriter(containedUnconfinedNetwork);
    // containedUnconfinedWriter.write("contained_unconfined", 2.5, 3.0);
    //
    // GmshWriter uncontainedConfinedWriter(uncontainedConfinedNetwork);
    // uncontainedConfinedWriter.write("uncontained_confined", 2.5);
    //
    // GmshWriter uncontainedUnconfinedWriter(uncontainedUnconfinedNetwork);
    // uncontainedUnconfinedWriter.write("uncontained_unconfined", 2.5);

    return 0;
}
