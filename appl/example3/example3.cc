#include <iostream>
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

    // we use the default sampler types -> normal distributions for all parameters
    using Distro = std::normal_distribution<ctype>;

    // Bounding box of the domain in which we want to place the entities
    const auto domainBBox = OCCUtilities::getBoundingBox(networkDomain);

    // sampler for disks (orientation 1)
    DiskSampler diskSampler(makeUniformPointSampler(domainBBox),                               // sampler for disk center points
                            Distro(30.0, 6.5),                        // major axis length: mean value & standard deviation
                            Distro(24.0, 4.5),                        // minor axis length: mean value & standard deviation
                            Distro(toRadians(0.0), toRadians(7.5)),   // rotation around x-axis: mean value & standard deviation
                            Distro(toRadians(0.0), toRadians(7.5)),   // rotation around y-axis: mean value & standard deviation
                            Distro(toRadians(0.0), toRadians(7.5)));  // rotation around z-axis: mean value & standard deviation

    // sampler for quadrilaterals (orientation 2)
    QuadrilateralSampler<3> quadSampler(makeUniformPointSampler(domainBBox),                               // sampler for quadrilateral center points
                                        Distro(toRadians(45.0), toRadians(5.0)),  // strike angle: mean value & standard deviation
                                        Distro(toRadians(90.0), toRadians(5.0)),  // dip angle: mean value & standard deviation
                                        Distro(45.0, 5.0),                        // edge length: mean value & standard deviation
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




    ///////////////////////////////////////////////////////////////////////
    // 3. Define constraints that should be fulfilled among the entities //
    //    of different orientations.                                     //
    ///////////////////////////////////////////////////////////////////////

    // Define constraints between entities of orientation 1
    using Constraints = EntityNetworkConstraints<ctype>;
    Constraints constraints1;
    constraints1.setMinDistance(2.5);
    constraints1.setMinIntersectingAngle(toRadians(25.0));
    constraints1.setMinIntersectionMagnitude(5.0);
    constraints1.setMinIntersectionDistance(2.5);

    // Define constraints between entities of orientation 2
    // we want to enforce larger spacing between those entities
    auto constraints2 = constraints1;
    constraints2.setMinDistance(5.0);

    // Define constraints between entities of different sets
    auto constraintsOnOther = constraints1;
    constraintsOnOther.setMinDistance(2.5);
    constraintsOnOther.setMinIntersectingAngle(toRadians(40.0));

    // We can use the constraints matrix to facilitate constraint evaluation
    EntityNetworkConstraintsMatrix<Constraints> constraintsMatrix;
    constraintsMatrix.addConstraints(constraints1,                  // constraint instance
                                     IdPair(diskSetId, diskSetId)); // set between which to use these constraints

    constraintsMatrix.addConstraints(constraints2,                  // constraint instance
                                     IdPair(quadSetId, quadSetId)); // set between which to use these constraints

    constraintsMatrix.addConstraints(constraintsOnOther,              // constraint instance
                                     {IdPair(diskSetId, quadSetId),   // sets between which to use these constraints
                                      IdPair(quadSetId, diskSetId)}); // sets between which to use these constraints

    // Moreover, we define constraints w.r.t. the domain boundary
    auto constraintsOnDomain = constraints1;
    constraintsOnDomain.setMinIntersectingAngle(toRadians(15.0));




    ///////////////////////////
    // 4. Network generation //
    ///////////////////////////

    // Helper class to store an arbitrary number
    // of entity sets of different geometry types.
    MultiGeometryEntitySet<Disk, Quad> entitySets;

    // Helper class for terminal output of the creation
    // progress and definition of stop criterion etc
    SamplingStatus status;
    status.setTargetCount(diskSetId, 12); // we want 11 entities of orientation 1
    status.setTargetCount(quadSetId, 16); // we want 13 entities of orientation 2

    // The actual network generation loop
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

        // Moreover, we want to avoid small fragments (< 250 m²)
        const auto containedArea = computeContainedMagnitude(geom, networkDomain);
        if (containedArea < 350.0)
        { status.increaseRejectedCounter(); continue; }

        // enforce constraints w.r.t. to the other entities
        if (!constraintsMatrix.evaluate(entitySets, geom, id))
        { status.increaseRejectedCounter(); continue; }

        // enforce constraints w.r.t. the domain boundaries
        if (!constraintsOnDomain.evaluate(domainBoundaryFaces, geom))
        { status.increaseRejectedCounter(); continue; }

        // the geometry is admissible
        entitySets.addEntity(geom, id);
        status.increaseCounter(id);
        status.print();

        // keep track of entity area
        containedNetworkArea += containedArea;
    }

    // print the final entity density
    const auto density = containedNetworkArea/domainVolume;
    std::cout << "\nEntity density of the contained network: " << density << " m²/m³" << std::endl;




    //////////////////////////////////////////////////////////////////////////
    // 5. The entities of the network have been created. We can now         //
    //    construct different types of networks (contained, confined, etc.) //
    //////////////////////////////////////////////////////////////////////////

    // construct and write a contained network, i.e. write out both network and domain.
    std::cout << "Building and writing contained, confined network" << std::endl;
    ContainedEntityNetworkBuilder<ctype> builder;

    // add sub-domains
    builder.addConfiningSubDomain(solids[0],     Id(1));
    builder.addConfiningSubDomain(networkDomain, Id(2));
    builder.addConfiningSubDomain(solids[2],     Id(3));

    // The entites that we created all belong to subdomain 2
    // Use the convenience function to add all entities to the network builder.
    // The last argument provides the id of the sub-domain to which the entities belong.
    entitySets.exportEntitySets(builder, Id(2));

    // Note that one can also call:
    // entitySets.exportEntitySets({diskSetId, quadSetId}, containedConfinedBuilder, Id(2));
    // This overload can be used to define specific entity sets that should be added to the builder

    // now we can build and write out the network in Gmsh file format
    GmshWriter gmshWriter(builder.build());
    gmshWriter.write("contained_confined", // body of the filename to be used (will add .geo)
                     2.5,                  // mesh size to be used on entities
                     5.0);                 // mesh size to be used on domain boundaries

    // we can also not confine the network to its sub-domain,
    // simply by adding the sub-domains as non-confining
    std::cout << "Building and writing contained, unconfined network" << std::endl;
    builder.clear();
    builder.addSubDomain(solids[0],     Id(1));
    builder.addSubDomain(networkDomain, Id(2));
    builder.addSubDomain(solids[2],     Id(3));
    entitySets.exportEntitySets(builder, Id(2));

    gmshWriter = GmshWriter(builder.build());
    gmshWriter.write("contained_unconfined", 2.5, 5.0);

    // We could also only write out the network, without the domain
    // For example, confining the network to the sub-domain...
    std::cout << "Building and writing uncontained, confined network" << std::endl;
    EntityNetworkBuilder<ctype> uncontainedBuilder;
    uncontainedBuilder.addConfiningSubDomain(networkDomain, Id(2));
    entitySets.exportEntitySets(uncontainedBuilder, Id(2));
    gmshWriter = GmshWriter(uncontainedBuilder.build());
    gmshWriter.write("uncontained_confined", 2.5);

    // ... or not confining it
    std::cout << "Building and writing uncontained, unconfined network" << std::endl;
    uncontainedBuilder.clear();
    entitySets.exportEntitySets(uncontainedBuilder);
    gmshWriter = GmshWriter(uncontainedBuilder.build());
    gmshWriter.write("uncontained_unconfined", 2.5);

    return 0;
}
