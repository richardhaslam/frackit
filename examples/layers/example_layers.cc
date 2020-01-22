#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <random>

// Contains functionality to read in shapes
#include <BRepTools.hxx>

// OpenCascade shape class
#include <TopoDS_Shape.hxx>

// basic geometry types
#include <frackit/geometry/box.hh>
#include <frackit/geometry/disk.hh>

// headers containing functions to compute lengths/areas/volumes
#include <frackit/magnitude/magnitude.hh>
#include <frackit/magnitude/containedmagnitude.hh>

// parser functions between internal & brep data structures
#include <frackit/occ/breputilities.hh>

// sampler for points and disks
#include <frackit/sampling/makeuniformpointsampler.hh>
#include <frackit/sampling/disksampler.hh>

// constraints to be enforced on the network (distance, angles, etc.)
#include <frackit/entitynetwork/constraints.hh>

// builder class for creating networks of entities confined in (sub-)domains
#include <frackit/entitynetwork/networkbuilder.hh>

// writes an entity network to a meshable Gmsh .geo file format
#include <frackit/io/gmshwriter.hh>

int main(int argc, char** argv)
{
    //! print welcome message
    std::cout << "\n\n"
              << "###################################################################\n"
              << "## Example: Creation of disk-shaped entities in a layered medium ##\n"
              << "###################################################################\n"
              << "\n\n";

    /////////////////////////////////////////////////////
    // 1. Read in the domain geometry from .brep file. //
    //    The file name is defined in CMakeLists.txt   //
    /////////////////////////////////////////////////////
    using namespace Frackit;
    const auto domainFileName = BREPFILE;
    const auto domainShape = OCCUtilities::readShape(domainFileName);

    // obtain the three solids contained in the file
    const auto solids = OCCUtilities::getSolids(domainShape);
    if (solids.size() != 3)
        throw std::runtime_error("Unexpected result read from .brep file");

    // the domain we want to create a network in is the center one
    const auto& networkDomain = solids[1];

    // compute domain volume and get its boundaries
    const auto domainVolume = computeMagnitude(networkDomain);
    const auto domainShells = OCCUtilities::getShells(networkDomain);

    if (domainShells.size() > 1)
        throw std::runtime_error(std::string("More than one shell read from solid"));
    else if (domainShells.empty())
        throw std::runtime_error(std::string("Could not read shell from solid"));

    // get the boundary faces of our domain of intererst
    const auto domainShellFaces = OCCUtilities::getFaces(domainShells[0]);

    //////////////////////////////////////////////////////////////////////
    // 2. Make sampler classes to randomly sample points (disk centers) //
    //    and disks (using the sampled center poitns)                   //
    //////////////////////////////////////////////////////////////////////
    using ctype = double;
    using Disk = Disk<ctype>;

    // sample points within bounding box of domain
    const auto domainBBox = OCCUtilities::getBoundingBox(networkDomain);

    // use default sampler: samples uniformly along each coordinate direction in the box
    auto pointSampler = makeUniformPointSampler(domainBBox);

    // mean major/minor axis length to be used in the disk samples
    const ctype meanMajAxisLength = 20.0;
    const ctype meanMinAxisLength = 15.0;

    // sampler for disks of orientation 1
    DiskSampler diskSampler1(pointSampler,                                                      // sampler for disk center points
                             std::normal_distribution<ctype>(meanMajAxisLength, 6.5),           // major axis length: mean value & standard deviation
                             std::normal_distribution<ctype>(meanMinAxisLength, 4.5),           // minor axis length: mean value & standard deviation
                             std::normal_distribution<ctype>(toRadians(0.0), toRadians(7.5)),   // rotation around x-axis: mean value & standard deviation
                             std::normal_distribution<ctype>(toRadians(0.0), toRadians(7.5)),   // rotation around y-axis: mean value & standard deviation
                             std::normal_distribution<ctype>(toRadians(0.0), toRadians(7.5)));  // rotation around z-axis: mean value & standard deviation

    // sampler for disks of orientation 2
    DiskSampler diskSampler2(pointSampler,                                                      // sampler for disk center points
                             std::normal_distribution<ctype>(meanMajAxisLength, 6.5),           // major axis length: mean value & standard deviation
                             std::normal_distribution<ctype>(meanMinAxisLength, 4.5),           // minor axis length: mean value & standard deviation
                             std::normal_distribution<ctype>(toRadians(90.0), toRadians(7.5)),  // rotation around x-axis: mean value & standard deviation
                             std::normal_distribution<ctype>(toRadians(0.0),  toRadians(7.5)),  // rotation around y-axis: mean value & standard deviation
                             std::normal_distribution<ctype>(toRadians(0.0),  toRadians(7.5))); // rotation around z-axis: mean value & standard deviation

    //! containers to store created entities
    std::vector<Disk> diskSet1;
    std::vector<Disk> diskSet2;

    //! enforce some constraints on the network
    auto constraintsOnSelf = makeDefaultConstraints<ctype>();
    auto constraintsOnOther = makeDefaultConstraints<ctype>();
    auto constraintsOnDomain = makeDefaultConstraints<ctype>();

    // constraints among disks of the same set
    constraintsOnSelf.setMinDistance(2.5);
    constraintsOnSelf.setMinIntersectingAngle(toRadians(45.0));
    constraintsOnSelf.setMinIntersectionMagnitude(2.5);
    constraintsOnSelf.setMinIntersectionDistance(2.0);

    // constraints with the other disk set
    constraintsOnOther.setMinDistance(1.5);
    constraintsOnOther.setMinIntersectingAngle(toRadians(10.0));
    constraintsOnOther.setMinIntersectionMagnitude(2.5);
    constraintsOnOther.setMinIntersectionDistance(2.0);

    // constraints with the domain boundary
    constraintsOnDomain.setMinDistance(1.5);
    constraintsOnDomain.setMinIntersectingAngle(toRadians(5.0));
    constraintsOnDomain.setMinIntersectionMagnitude(20.0);
    constraintsOnDomain.setMinIntersectionDistance(2.0);

    //! create pseudo-random number generator to select set during creation
    std::default_random_engine randomNumber{std::random_device{}()};

    //! keep track of number of created & accepted entities
    std::size_t total = 0;
    std::size_t accepted = 0;
    std::size_t accepted_1 = 0;
    std::size_t accepted_2 = 0;

    //! Define how many entities should be created
    const std::size_t numTargetEntities_1 = 4;
    const std::size_t numTargetEntities_2 = 4;

    //! print header of status output
    std::cout << "##############################################################################################################################\n"
              << "# No. of entities   |   No. rejected entities   |   Acceptance Ratio   |   Current entity density [m²/m³]   |   Progress [%] #\n"
              << "#----------------------------------------------------------------------------------------------------------------------------#\n";

    ctype currentDensity = 0.0;
    ctype currentDiskArea = 0.0;
    while (accepted_1 != numTargetEntities_1 || accepted_2 != numTargetEntities_2)
    {
        bool createSecondary = randomNumber()%2;
        if (!createSecondary && accepted_1 == numTargetEntities_1)
            createSecondary = true;
        else if (createSecondary && accepted_2 == numTargetEntities_2)
            createSecondary = false;

        auto disk = createSecondary ? diskSampler1()
                                    : diskSampler2();
        total++;

        // We don't want ellipses of too large aspect ratio
        if (disk.majorAxisLength() > disk.minorAxisLength()*3.0) continue;

        auto& diskSetSelf = createSecondary ? diskSet2 : diskSet1;
        const auto& diskSetOther = createSecondary ? diskSet1 : diskSet2;

        // enforce constraints w.r.t. to the own disk set
        if (!constraintsOnSelf.evaluate(diskSetSelf, disk)) continue;
        // enforce constraints w.r.t. to the other disk set
        if (!constraintsOnOther.evaluate(diskSetOther, disk)) continue;
        // enforce constraints w.r.t. the domain boundaries
        if (!constraintsOnDomain.evaluate(domainShellFaces, disk)) continue;

        // reject if disk area contained in the domain is too small
        const auto containedArea = computeContainedMagnitude(disk, networkDomain);
        if (containedArea < meanMajAxisLength*meanMinAxisLength*M_PI/5.0) continue;

        // if we get here, the disk is admissible
        diskSetSelf.push_back(std::move(disk));
        accepted++;

        if (createSecondary) accepted_2++;
        else accepted_1++;

        // compute new density (use minimum w.r.t. domain/network volume)
        currentDiskArea += containedArea;
        currentDensity = currentDiskArea/domainVolume;
        std::cout << "  "   << std::setw(18) << std::setfill(' ')
                            << std::to_string(diskSet1.size() + diskSet2.size()) + std::string(9, ' ')
                  << "|   " << std::setprecision(6) << std::setw(24) << std::setfill(' ')
                            << std::to_string(total-accepted) + std::string(12, ' ')
                  << "|   " << std::setprecision(6) << std::setw(19) << std::setfill(' ')
                            << std::to_string( double(double(accepted)/double(total)) ) + std::string(7, ' ')
                  << "|   " << std::setprecision(6) << std::setw(33) << std::setfill(' ')
                            << std::to_string(currentDensity) + std::string(17, ' ')
                  << "|   " << std::string(3, ' ') + std::to_string( double(accepted_1+accepted_2)/double(numTargetEntities_1+numTargetEntities_2) ) << std::endl;
    }


    //////////////////////////////////////////////////////////////////////////
    // 3. The entities of the network have been created. We can now         //
    //    construct different types of networks (contained, confined, etc.) //
    //////////////////////////////////////////////////////////////////////////

    std::cout << "Constructing contained, confined network" << std::endl;
    // 3.1 Construct a network such that the entities are confined in the
    //     sub-domain and which contains information about the containing domain.
    ContainedEntityNetworkBuilder containedConfinedBuilder;

    // define sub-domains
    containedConfinedBuilder.addConfiningSubDomain(solids[0],     Id(1));
    containedConfinedBuilder.addConfiningSubDomain(networkDomain, Id(2));
    containedConfinedBuilder.addConfiningSubDomain(solids[2],     Id(3));

    // define entity network for sub-domain 2
    containedConfinedBuilder.addSubDomainEntities(diskSet1, Id(2));
    containedConfinedBuilder.addSubDomainEntities(diskSet2, Id(2));

    // build network
    const auto containedConfinedNetwork = containedConfinedBuilder.build();

    std::cout << "Constructing contained, unconfined network" << std::endl;
    // 3.2 Construct a network such that the entities are NOT confined in the
    //     sub-domain and which contains information about the containing domain.
    ContainedEntityNetworkBuilder containedUnconfinedBuilder;

    // define sub-domains
    containedUnconfinedBuilder.addSubDomain(solids[0],     Id(1));
    containedUnconfinedBuilder.addSubDomain(solids[2],     Id(3));
    containedUnconfinedBuilder.addSubDomain(networkDomain, Id(2));

    // define entity network for sub-domain 2
    containedUnconfinedBuilder.addSubDomainEntities(diskSet1, Id(2));
    containedUnconfinedBuilder.addSubDomainEntities(diskSet2, Id(2));

    // build network
    const auto containedUnconfinedNetwork = containedUnconfinedBuilder.build();

    std::cout << "Constructing uncontained, confined network" << std::endl;
    // 3.3 Construct a network such that the entities are confined in the
    //     sub-domain, but, which does not contain information about the domains.
    //     This is computationally more efficient if only the 2d network is of interest.
    //     Moreover, the files written to disk are smaller as the domain is not contained.
    EntityNetworkBuilder uncontainedConfinedBuilder;

    // define sub-domains
    uncontainedConfinedBuilder.addConfiningSubDomain(networkDomain, Id(1));

    // define entity network for sub-domain 2
    uncontainedConfinedBuilder.addSubDomainEntities(diskSet1, Id(1));
    uncontainedConfinedBuilder.addSubDomainEntities(diskSet2, Id(1));

    // build network
    const auto uncontainedConfinedNetwork = uncontainedConfinedBuilder.build();

    std::cout << "Constructing uncontained, unconfined network" << std::endl;
    // 3.4 Construct the network only, independent of any domains
    EntityNetworkBuilder uncontainedUnconfinedBuilder;
    uncontainedUnconfinedBuilder.addEntities(diskSet1);
    uncontainedUnconfinedBuilder.addEntities(diskSet2);

    // build network
    const auto uncontainedUnconfinedNetwork = uncontainedUnconfinedBuilder.build();


    ///////////////////////////////////////////////////////
    // 4. Write the networks to disk in Gmsh .geo format //
    ///////////////////////////////////////////////////////
    std::cout << "Write networks to disk" << std::endl;
    GmshWriter containedConfinedWriter(containedConfinedNetwork);
    containedConfinedWriter.write("contained_confined", // body of the filename to be used (will add .geo)
                                  2.5,                  // mesh size to be used on entities
                                  5.0);                 // mesh size to be used on domain boundaries

    GmshWriter containedUnconfinedWriter(containedUnconfinedNetwork);
    containedUnconfinedWriter.write("contained_unconfined", 2.5, 3.0);

    GmshWriter uncontainedConfinedWriter(uncontainedConfinedNetwork);
    uncontainedConfinedWriter.write("uncontained_confined", 2.5);

    GmshWriter uncontainedUnconfinedWriter(uncontainedUnconfinedNetwork);
    uncontainedUnconfinedWriter.write("uncontained_unconfined", 2.5);

    return 0;
}
