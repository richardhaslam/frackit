#include <iostream>
#include <random>

#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Compound.hxx>

#include <frackit/geometry/cylinder.hh>
#include <frackit/precision/defaultepsilon.hh>
#include <frackit/magnitude/containedmagnitude.hh>
#include <frackit/occ/breputilities.hh>

#include <frackit/sampling/geometrypointsampler.hh>
#include <frackit/sampling/disksampler.hh>

#include <frackit/entitynetwork/constraints.hh>

//! create a disk network embedded in a cylinder
int main()
{
    //! print welcome message
    std::cout << "\n\n"
              << "##########################################\n"
              << "## Creating a confined network of disks ##\n"
              << "##########################################\n"
              << "\n\n";

    using namespace Frackit;
    using ctype = double;
    using Disk = Disk<ctype>;
    using Domain = Cylinder<ctype>;

    Domain domain(0.5, 1.0);
    auto pointSampler = makeUniformPointSampler(domain);

    // sampler for disks of orientation 1
    DiskSampler diskSampler_1(std::normal_distribution<ctype>(0.35, 0.1),
                              std::normal_distribution<ctype>(0.225, 0.05),
                              std::normal_distribution<ctype>(toRadians(25.0), toRadians(5.0)),
                              std::normal_distribution<ctype>(toRadians(0.0),  toRadians(5.0)),
                              std::normal_distribution<ctype>(toRadians(45.0), toRadians(5.0)));

    // sampler for disks of orientation 1
    DiskSampler diskSampler_2(std::normal_distribution<ctype>(0.35, 0.1),
                              std::normal_distribution<ctype>(0.225, 0.05),
                              std::normal_distribution<ctype>(toRadians(-35.0), toRadians(5.0)),
                              std::normal_distribution<ctype>(toRadians(0.0),   toRadians(5.0)),
                              std::normal_distribution<ctype>(toRadians(-45.0), toRadians(5.0)));

    //! containers to store created entities
    std::vector<Disk> diskSet1;
    std::vector<Disk> diskSet2;

    //! enforce some constraints on the network
    auto constraintsOnSelf = makeDefaultConstraints<ctype>();
    auto constraintsOnOther = makeDefaultConstraints<ctype>();
    auto constraintsOnDomain = makeDefaultConstraints<ctype>();

    // constraints among disks of the same set
    constraintsOnSelf.setMinDistance(0.1);
    constraintsOnSelf.setMinIntersectingAngle(toRadians(45.0));
    constraintsOnSelf.setMinIntersectionMagnitude(0.01);
    constraintsOnSelf.setMinIntersectionDistance(0.05);

    // constraints with the other disk set
    constraintsOnOther.setMinDistance(0.01);
    constraintsOnOther.setMinIntersectingAngle(toRadians(15.0));
    constraintsOnOther.setMinIntersectionMagnitude(0.01);
    constraintsOnOther.setMinIntersectionDistance(0.05);

    // constraints with the other disk set
    constraintsOnDomain.setMinDistance(0.01);
    constraintsOnDomain.setMinIntersectingAngle(toRadians(15.0));
    constraintsOnDomain.setMinIntersectionMagnitude(0.01);
    constraintsOnDomain.setMinIntersectionDistance(0.05);

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
    const std::size_t numTargetEntities_1 = 5;
    const std::size_t numTargetEntities_2 = 5;
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
        if (!constraintsOnDomain.evaluate(domain.bottomFace(), disk)) continue;
        // enforce constraints w.r.t. the domain boundaries
        if (!constraintsOnDomain.evaluate(domain.topFace(), disk)) continue;
        // enforce constraints w.r.t. the domain boundaries
        if (!constraintsOnDomain.evaluate(domain.lateralFace(), disk)) continue;

        // reject if intersection with domain is too small (here: 0.2m²)
        const auto containedArea = computeContainedMagnitude(disk, domain);
        if (containedArea < 0.1) continue;

        // if we get here, the disk is admissible
        diskSetSelf.push_back(std::move(disk));
        accepted++;

        if (createSecondary) accepted_2++;
        else accepted_1++;

        // compute new density (use minimum w.r.t. domain/network volume)
        currentDiskArea += containedArea;
        currentDensity = currentDiskArea/domain.volume();
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
    b.Add(c, OCCUtilities::getShape(domain));
    BRepTools::Write(c, "disknetwork.brep");

    // fragment the network and write in separate file
    std::vector<TopoDS_Shape> allDisks;
    for (const auto& disk : diskSet1) allDisks.push_back(OCCUtilities::getShape(disk));
    for (const auto& disk : diskSet2) allDisks.push_back(OCCUtilities::getShape(disk));

    const auto fragmentedNetwork = OCCUtilities::fragment(allDisks, defaultEpsilon(domain));
    BRepTools::Write(fragmentedNetwork, "disknetwork_fragmented.brep");

    const auto confinedNetwork = OCCUtilities::intersect(fragmentedNetwork,
                                                         OCCUtilities::getShape(domain),
                                                         defaultEpsilon(domain));
    BRepTools::Write(confinedNetwork, "disknetwork_confined.brep");

    std::vector<TopoDS_Shape> allShapes({confinedNetwork, OCCUtilities::getShape(domain)});
    BRepTools::Write(OCCUtilities::fragment(allShapes, defaultEpsilon(domain)), "disknetwork_medium.brep");

    return 0;
}
