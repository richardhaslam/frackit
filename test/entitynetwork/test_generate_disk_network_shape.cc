#include <iostream>
#include <random>
#include <fstream>

#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Solid.hxx>

#include <frackit/geometry/box.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/precision/precision.hh>
#include <frackit/magnitude/magnitude.hh>
#include <frackit/magnitude/containedmagnitude.hh>
#include <frackit/occ/breputilities.hh>

#include <frackit/sampling/pointsampler.hh>
#include <frackit/sampling/disksampler.hh>

#include <frackit/entitynetwork/constraints.hh>

//! create a disk network embedded in a solid shape
int main(int argc, char** argv)
{
    //! print welcome message
    std::cout << "\n\n"
              << "##########################################\n"
              << "## Creating a confined network of disks ##\n"
              << "##########################################\n"
              << "\n\n";

    // brep file name is specified in CMakeLists.txt
    std::ifstream domainShapeFile(BREPFILE);
    if (!domainShapeFile)
        throw std::runtime_error("Could not open shape file");

    TopoDS_Shape domainShape;
    BRepTools::Read(domainShape, domainShapeFile, BRep_Builder());

    // obtain the (single) solid contained in the file & get its boundaries
    using namespace Frackit;
    const auto solids = OCCUtilities::getSolids(domainShape);

    assert(solids.size() == 1);
    const auto& domain = solids[0];
    const auto domainVolume = computeMagnitude(domain);
    const auto domainShells = OCCUtilities::getShells(domain);
    assert(domainShells.size() == 1);
    const auto domainShellFaces = OCCUtilities::getFaces(domainShells[0]);

    // create the disk samplers
    using ctype = double;
    using Disk = Disk<ctype>;

    // sample points within bounding box of domain
    auto pointSampler = makeUniformPointSampler(OCCUtilities::getBoundingBox(domain));

    // sampler for disks of orientation 1
    DiskSampler diskSampler_1(std::normal_distribution<ctype>(40.0, 6.5),
                              std::normal_distribution<ctype>(20.0, 4.5),
                              std::normal_distribution<ctype>(toRadians(0.0), toRadians(5.0)),
                              std::normal_distribution<ctype>(toRadians(0.0), toRadians(5.0)),
                              std::normal_distribution<ctype>(toRadians(0.0), toRadians(5.0)));

    // sampler for disks of orientation 1
    DiskSampler diskSampler_2(std::normal_distribution<ctype>(40.0, 6.5),
                              std::normal_distribution<ctype>(20.0, 4.5),
                              std::normal_distribution<ctype>(toRadians(45.0), toRadians(5.0)),
                              std::normal_distribution<ctype>(toRadians(0.0),  toRadians(5.0)),
                              std::normal_distribution<ctype>(toRadians(0.0),  toRadians(5.0)));

    //! containers to store created entities
    std::vector<Disk> diskSet1;
    std::vector<Disk> diskSet2;

    //! enforce some constraints on the network
    auto constraintsOnSelf = makeDefaultConstraints<ctype>();
    auto constraintsOnOther = makeDefaultConstraints<ctype>();
    auto constraintsOnDomain = makeDefaultConstraints<ctype>();

    // constraints among disks of the same set
    constraintsOnSelf.setMinDistance(10.0);
    constraintsOnSelf.setMinIntersectingAngle(toRadians(45.0));
    constraintsOnSelf.setMinIntersectionMagnitude(2.5);
    constraintsOnSelf.setMinIntersectionDistance(2.0);

    // constraints with the other disk set
    constraintsOnOther.setMinDistance(2.5);
    constraintsOnOther.setMinIntersectingAngle(toRadians(10.0));
    constraintsOnOther.setMinIntersectionMagnitude(10.0);
    constraintsOnOther.setMinIntersectionDistance(2.0);

    // constraints with the domain
    constraintsOnDomain.setMinDistance(2.5);
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

    //! print header of status output
    std::cout << "##############################################################################################################################\n"
              << "# No. of entities   |   No. rejected entities   |   Acceptance Ratio   |   Current entity density [m²/m³]   |   Progress [%] #\n"
              << "#----------------------------------------------------------------------------------------------------------------------------#\n";

    ctype currentDensity = 0.0;
    ctype currentDiskArea = 0.0;
    const std::size_t numTargetEntities_1 = 3;
    const std::size_t numTargetEntities_2 = 3;
    while (accepted_1 != numTargetEntities_1 || accepted_2 != numTargetEntities_2)
    {
        bool createSecondary = randomNumber()%2;
        if (!createSecondary && accepted_1 == numTargetEntities_1)
            createSecondary = true;
        else if (createSecondary && accepted_2 == numTargetEntities_2)
            createSecondary = false;

        auto disk = createSecondary ? diskSampler_1( pointSampler() )
                                    : diskSampler_2( pointSampler() );
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
        bool violates = false;
        for (const auto& face : domainShellFaces)
            if (!constraintsOnDomain.evaluate(face, disk))
            { violates = true; break; }
        if (violates) continue;

        // reject if intersection with domain is too small (here: 0.2m²)
        const auto containedArea = computeContainedMagnitude(disk, domain);
        if (containedArea < 20.0*10.0*M_PI/5.0) continue;

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

    // create a single compound and write to .brep file
    BRep_Builder b;
    TopoDS_Compound c;
    b.MakeCompound(c);
    for (const auto& disk : diskSet1) b.Add(c, OCCUtilities::getShape(disk));
    for (const auto& disk : diskSet2) b.Add(c, OCCUtilities::getShape(disk));
    b.Add(c, domain);
    BRepTools::Write(c, "disknetwork.brep");

    // fragment the network and write in separate file
    std::vector<TopoDS_Shape> allDisks;
    for (const auto& disk : diskSet1) allDisks.push_back(OCCUtilities::getShape(disk));
    for (const auto& disk : diskSet2) allDisks.push_back(OCCUtilities::getShape(disk));

    const auto fragmentedNetwork = OCCUtilities::fragment(allDisks, Frackit::Precision<ctype>::confusion());
    BRepTools::Write(fragmentedNetwork, "disknetwork_fragmented.brep");

    const auto confinedNetwork = OCCUtilities::intersect(fragmentedNetwork,
                                                         domain,
                                                         Frackit::Precision<ctype>::confusion());
    BRepTools::Write(confinedNetwork, "disknetwork_confined.brep");

    std::vector<TopoDS_Shape> allShapes({confinedNetwork, domain});
    BRepTools::Write(OCCUtilities::fragment(allShapes, Frackit::Precision<ctype>::confusion()),
                    "disknetwork_medium.brep");

    std::cout << "All tests passed" << std::endl;

    return 0;
}
