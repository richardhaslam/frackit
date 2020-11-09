#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <random>
#include <vector>
#include <cmath>
#include <chrono>

// utility functions
#include <frackit/common/math.hh>

// basic geometry types
#include <frackit/geometry/box.hh>
#include <frackit/geometry/sphere.hh>
#include <frackit/geometry/polygon.hh>
#include <frackit/geometry/cylinder.hh>
#include <frackit/geometry/hollowcylinder.hh>
#include <frackit/geometryutilities/getboundingbox.hh>

// headers containing functions to compute lengths/areas/volumes/distances
#include <frackit/distance/distance.hh>
#include <frackit/magnitude/magnitude.hh>

// sampler for points and polygons
#include <frackit/sampling/makeuniformpointsampler.hh>
#include <frackit/sampling/polygonsampler.hh>
#include <frackit/sampling/status.hh>

// constraints to be enforced on the network (distance, angles, etc.)
#include <frackit/entitynetwork/constraints.hh>

// brep utilities of OpenCascade, for boolean operations on shapes
#include <TopoDS_Face.hxx>
#include <frackit/occ/breputilities.hh>

// intersection detection for internal geometry types
#include <frackit/intersection/intersect.hh>
#include <frackit/intersection/emptyintersection.hh>

// builder class for creating networks of entities confined in (sub-)domains
#include <frackit/entitynetwork/networkbuilder.hh>

// writes an entity network to a meshable Gmsh .geo file format
#include <frackit/io/gmshwriter.hh>

int main(int argc, char** argv)
{
    // get start time to compute CPU time later
    const auto start = std::chrono::steady_clock::now();

    //! print welcome message
    std::cout << "\n\n"
              << "#########################################################################\n"
              << "## Example 4: Creation of polygonal entities around a tunnel structure ##\n"
              << "#########################################################################\n"
              << "\n\n\n";

    // Define some types used here
    using namespace Frackit;
    using ctype = double;
    using Box = Frackit::Box<ctype>;
    using Sphere = Frackit::Sphere<ctype>;
    using Cylinder = Frackit::Cylinder<ctype>;
    using HollowCylinder = Frackit::HollowCylinder<ctype>;
    using Polygon = Frackit::Polygon<ctype, 3>;

    using Point = typename Polygon::Point;
    using Circle = typename Cylinder::Circle;
    using Direction = typename Cylinder::Direction;
    using Vector = typename Direction::Vector;

    /////////////////////////////
    // 1. Construct the domain //
    /////////////////////////////

    // In this example, we consider a fraction of a geological layer
    // around a tunnel structure. To this end, we construct the geological
    // layer as box, and intersect it with a sphere in order to restrict
    // the model domain to a smaller region.

    // the layer has a thickness of 40m, and initially we consider 100m lateral extent
    const Box layer(-50, -50, -20.0, 50, 50, 20.0);
    const Sphere domainSphere(Point(0.0, 0.0, 0.0), /*radius*/25.0);

    // a tunnel with a radius of 5m crosses the domain
    const ctype tunnelRadius = 5.0;
    const Direction tunnelDirection(Vector(1.0, 0.0, 0.0));
    const Circle tunnelBase(Point(-25.0, 0.0, 0.0), tunnelDirection, tunnelRadius);
    const Cylinder tunnel(tunnelBase, /*height*/50.0);

    // the domain is the layer constrained to the domainSphere and cut by the tunnel
    const auto layerSphereCut = OCCUtilities::intersect(OCCUtilities::getShape(layer),
                                                        OCCUtilities::getShape(domainSphere),
                                                        /*eps*/1e-6);
    auto solids = OCCUtilities::getSolids(layerSphereCut);
    if (solids.size() != 1) throw std::runtime_error("Intersection operation with domainSphere should yield a single solid");

    const auto tunnelCut = OCCUtilities::cut(solids[0], OCCUtilities::getShape(tunnel), /*eps*/1e-6);
    solids = OCCUtilities::getSolids(tunnelCut);
    if (solids.size() != 1) throw std::runtime_error("Cut operation with tunnel should yield a single solid");

    // get the domain shape
    const auto& domain = solids[0];

    // get the OCC reprensentation of the sphere surface
    const auto sphereFaces = OCCUtilities::getFaces(OCCUtilities::getShape(domainSphere));
    if (sphereFaces.size() != 1) throw std::runtime_error("We expect the sphere to consist of a single face");
    const auto& sphereSurfaceShape = sphereFaces[0];

    // get the OCC representation of tunnel surface
    const auto tunnelSurfaceShape = OCCUtilities::getShape(tunnel.lateralFace());


    ////////////////////////////////////////////////////////////////////////
    // 2. Make sampler classes to randomly sample points (entity centers) //
    //    and polygons as entities (using the sampled center points)      //
    ////////////////////////////////////////////////////////////////////////

    // we use the default sampler types, thus, default distributions (see traits classes)
    // This means, normal distributions for angles and uniform distributions for the size
    using PolygonSampler = Frackit::PolygonSampler<3, ctype>;

    // lambda function to construct polygon sampler with given point sampler & orientation
    auto getSampler = [] (const auto& pointSampler, unsigned int orientation) -> PolygonSampler
    {
        using NormalDistro = std::normal_distribution<ctype>;
        auto sizeDistro = std::uniform_real_distribution<ctype>(3.0, 6.0);
        auto numCornerDistro = std::uniform_int_distribution<int>(4, 9);

        switch (orientation)
        {
            case 1:
                return PolygonSampler(pointSampler,                                   // sampler for polygon centers
                                      NormalDistro(toRadians(15.0), toRadians(10.0)), // strike angle: mean value & standard deviation
                                      NormalDistro(toRadians(90.0), toRadians(10.0)), // dip angle: mean value & standard deviation
                                      /*strikeLength/*/sizeDistro, /*dipLength*/sizeDistro, numCornerDistro);
            case 2:
                return PolygonSampler(pointSampler,                                  // sampler for polygon centers
                                      NormalDistro(toRadians(0.0), toRadians(10.0)), // strike angle: mean value & standard deviation
                                      NormalDistro(toRadians(0.0), toRadians(10.0)), // dip angle: mean value & standard deviation
                                      /*strikeLength/*/sizeDistro, /*dipLength*/sizeDistro, numCornerDistro);
            case 3:
                return PolygonSampler(pointSampler,                                    // sampler for polygon centers
                                      NormalDistro(toRadians(-75.0), toRadians(10.0)), // strike angle: mean value & standard deviation
                                      NormalDistro(toRadians(90.0), toRadians(10.0)),  // dip angle: mean value & standard deviation
                                      /*strikeLength/*/sizeDistro, /*dipLength*/sizeDistro, numCornerDistro);
            default:
                throw std::runtime_error("Unsupported orientation index");
        }
    };


    ///////////////////////////////////////////////////////////////////////
    // 3. Define constraints that should be fulfilled among the entities //
    //    of different orientations.                                     //
    ///////////////////////////////////////////////////////////////////////

    // Define constraints between entities of orientation 1
    using Constraints = EntityNetworkConstraints<ctype>;
    Constraints constraints1;
    constraints1.setMinDistance(1.0);
    constraints1.setMinIntersectingAngle(toRadians(25.0));
    constraints1.setMinIntersectionMagnitude(0.1);
    constraints1.setMinIntersectionDistance(0.1);

    // Define constraints between entities of orientation 2
    // we want to enforce larger spacing between those entities
    auto constraints2 = constraints1;
    constraints2.setMinDistance(0.5);

    // Define constraints between entities of orientation 3
    auto constraints3 = constraints1;
    constraints2.setMinDistance(0.15);

    // Define constraints between entities of different sets
    auto constraintsOnOther = constraints1;
    constraintsOnOther.setMinDistance(0.25);
    constraintsOnOther.setMinIntersectingAngle(toRadians(30.0));

    // Define constraints to the domain boundary
    auto constraintsOnBoundary = constraints1;
    constraintsOnBoundary.setMinIntersectingAngle(toRadians(10.0));

    // When polygons intersect, the intersection might be very close to
    // one of the corner points. We want to avoid small length scales caused
    // by this, and this lambda provides the check for it.
    auto producesSmallLengthScale = [] (const auto& geometry, const auto& is) -> bool
    {
        // lambda to check if the intersection is too close (< 1mm) to a vertex of the geometry
        auto isTooClose = [&is] (const auto& v) { return computeDistance(v, is) < 1e-3; };

        const auto shape = OCCUtilities::getShape(geometry);
        const auto vertices = OCCUtilities::getVertices(shape);
        return std::any_of(vertices.begin(), vertices.end(), isTooClose);
    };


    ///////////////////////////
    // 4. Network generation //
    ///////////////////////////
    const std::size_t numTargetEntities = 75;
    std::vector<TopoDS_Face> entitiesSet1; entitiesSet1.reserve(numTargetEntities);
    std::vector<TopoDS_Face> entitiesSet2; entitiesSet2.reserve(numTargetEntities);
    std::vector<TopoDS_Face> entitiesSet3; entitiesSet3.reserve(numTargetEntities);

    // Define ids for the two entity sets
    const Id idSet1(1);
    const Id idSet2(2);
    const Id idSet3(3);

    // In a first step, we want to populate entities of orientation 1 around
    // the tunnel, and require each entity to be intersecting with it. To this end,
    // we sample center points within a hollow cylinder around the tunnel.
    SamplingStatus status; status.setTargetCount(idSet1, numTargetEntities);
    const HollowCylinder sampleCylinder1(tunnel.bottomFace().center(),
                                         tunnel.direction(),
                                         tunnel.radius(),
                                         tunnel.radius() + 2.0,
                                         tunnel.height());

    ctype area = 0.0;
    auto polygonSampler1 = getSampler(makeUniformPointSampler(sampleCylinder1), /*orientationIdx*/1);
    while (!status.finished())
    {
        const auto candidate = polygonSampler1();
        const auto is = intersect(tunnel.lateralFace(), candidate);

        // we only want entities that intersect the tunnel
        if (isEmptyIntersection(is))
        { status.increaseRejectedCounter(); continue; }

        // avoid small length scales produced by the intersection
        if (producesSmallLengthScale(candidate, is))
        { status.increaseRejectedCounter(); continue; }

        // compute the part contained in the domain
        const auto containedShape = OCCUtilities::intersect(OCCUtilities::getShape(candidate), domain, /*eps*/1e-6);
        const auto containedFaces = OCCUtilities::getFaces(containedShape);

        // skip everything that does not lead to a single contained face
        if (containedFaces.size() != 1)
        { status.increaseRejectedCounter(); continue; }

        // if too little of the candidate is inside the domain, reject it
        auto containedFace = containedFaces[0];
        const auto containedFaceArea = computeMagnitude(containedFace);
        if (containedFaceArea < 0.25*candidate.area())
        { status.increaseRejectedCounter(); continue; }

        // check constraints to other entities and domain boundary
        if (!constraints1.evaluate(entitiesSet1, containedFace))
        { status.increaseRejectedCounter(); continue; }
        if (!constraintsOnBoundary.evaluate(tunnelSurfaceShape, containedFace))
        { status.increaseRejectedCounter(); continue; }
        if (!constraintsOnBoundary.evaluate(sphereSurfaceShape, containedFace))
        { status.increaseRejectedCounter(); continue; }

        area += containedFaceArea;
        entitiesSet1.emplace_back(std::move(containedFace));
        status.increaseCounter(idSet1);
        status.print();
    }

    std::cout << "\n\n\n"
              << "## Step 2: for each entity so far, try to generate an intersecting entity with orientation 2\n"
              << "\n\n";

    // after rejecting 500 candidates for an entity, we move to the next one
    const std::size_t maxTriesForEntity = 500;
    auto printSkipMessage = [&maxTriesForEntity] ()
    {
        std::cout << "\t\t -- Skipped search for this entity after "
                  << maxTriesForEntity << " unsuccessful tries. --" << std::endl;
    };

    status.reset();
    status.setTargetCount(idSet2, entitiesSet1.size());
    for (const auto& e : entitiesSet1)
    {
        using std::sqrt;

        // sample an entity of orientation 2 in the vicinity of e
        const auto rawBBox = getBoundingBox(e);
        const auto charLength = sqrt(computeMagnitude(e));
        const Box sampleBox(rawBBox.xMin()-charLength, rawBBox.yMin()-charLength, rawBBox.zMin()-charLength,
                            rawBBox.xMax()+charLength, rawBBox.yMax()+charLength, rawBBox.zMax()+charLength);

        bool accepted = false;
        std::size_t tryCount = 0;
        auto sampler = getSampler(makeUniformPointSampler(sampleBox), /*orientationIdx*/2);

        while (!accepted && tryCount <= maxTriesForEntity)
        {
            tryCount++;
            const auto candidate = sampler();
            const auto is = intersect(candidate, e);

            // we only want intersecting entities
            if (isEmptyIntersection(is))
            { status.increaseRejectedCounter(); continue; }

            // avoid small length scales produced by the intersection
            if (producesSmallLengthScale(candidate, is))
            { status.increaseRejectedCounter(); continue; }
            if (producesSmallLengthScale(e, is))
            { status.increaseRejectedCounter(); continue; }

            // compute the part contained in the domain
            const auto containedShape = OCCUtilities::intersect(OCCUtilities::getShape(candidate), domain, /*eps*/1e-6);
            const auto containedFaces = OCCUtilities::getFaces(containedShape);

            // skip everything that does not lead to a single contained face
            if (containedFaces.size() != 1)
            { status.increaseRejectedCounter(); continue; }

            // if too little of the candidate is inside the domain, reject it
            auto containedFace = containedFaces[0];
            const auto containedFaceArea = computeMagnitude(containedFace);
            if (containedFaceArea < 0.25*candidate.area())
            { status.increaseRejectedCounter(); continue; }

            // check constraints to other entities and domain boundary
            if (!constraintsOnOther.evaluate(entitiesSet1, containedFace))
            { status.increaseRejectedCounter(); continue; }
            if (!constraints2.evaluate(entitiesSet2, containedFace))
            { status.increaseRejectedCounter(); continue; }
            if (!constraintsOnBoundary.evaluate(sphereSurfaceShape, containedFace))
            { status.increaseRejectedCounter(); continue; }
            if (!constraintsOnBoundary.evaluate(tunnelSurfaceShape, containedFace))
            { status.increaseRejectedCounter(); continue; }

            accepted = true;
            area += containedFaceArea;
            entitiesSet2.emplace_back(std::move(containedFace));
            status.increaseCounter(idSet2);
            status.print();
        }

        if (tryCount >= maxTriesForEntity) printSkipMessage();
    }

    std::cout << "\n\n\n"
              << "## Step 3: for each entity of orientation 2, try to generate an intersecting entity with orientation 3\n"
              << "\n\n";

    status.reset();
    status.setTargetCount(idSet3, entitiesSet2.size());
    for (const auto& e : entitiesSet2)
    {
        using std::sqrt;

        // sample an entity of orientation 2 in the vicinity of e
        const auto rawBBox = getBoundingBox(e);
        const auto charLength = sqrt(computeMagnitude(e));
        const Box sampleBox(rawBBox.xMin()-charLength, rawBBox.yMin()-charLength, rawBBox.zMin()-charLength,
                            rawBBox.xMax()+charLength, rawBBox.yMax()+charLength, rawBBox.zMax()+charLength);

        auto sampler = getSampler(makeUniformPointSampler(sampleBox), /*orientationIdx*/3);
        bool accepted = false;
        std::size_t tryCount = 0;

        while (!accepted && tryCount <= maxTriesForEntity)
        {
            tryCount++;
            auto candidate = sampler();
            const auto is = intersect(candidate, e);

            // we only want intersecting entities
            if (isEmptyIntersection(is))
            { status.increaseRejectedCounter(); continue; }

            // avoid small length scales produced by the intersection
            if (producesSmallLengthScale(candidate, is))
            { status.increaseRejectedCounter(); continue; }
            if (producesSmallLengthScale(e, is))
            { status.increaseRejectedCounter(); continue; }

            // compute the part contained in the domain
            const auto containedShape = OCCUtilities::intersect(OCCUtilities::getShape(candidate), domain, /*eps*/1e-6);
            const auto containedFaces = OCCUtilities::getFaces(containedShape);

            // skip everything that does not lead to a single contained face
            if (containedFaces.size() != 1)
            { status.increaseRejectedCounter(); continue; }

            // if too little of the candidate is inside the domain, reject it
            auto containedFace = containedFaces[0];
            const auto containedFaceArea = computeMagnitude(containedFace);
            if (containedFaceArea < 0.25*candidate.area())
            { status.increaseRejectedCounter(); continue; }

            // check constraints to other entities and domain boundary
            if (!constraintsOnOther.evaluate(entitiesSet1, containedFace))
            { status.increaseRejectedCounter(); continue; }
            if (!constraintsOnOther.evaluate(entitiesSet2, containedFace))
            { status.increaseRejectedCounter(); continue; }
            if (!constraints3.evaluate(entitiesSet3, containedFace))
            { status.increaseRejectedCounter(); continue; }
            if (!constraintsOnBoundary.evaluate(sphereSurfaceShape, containedFace))
            { status.increaseRejectedCounter(); continue; }
            if (!constraintsOnBoundary.evaluate(tunnelSurfaceShape, containedFace))
            { status.increaseRejectedCounter(); continue; }

            accepted = true;
            area += containedFaceArea;
            entitiesSet3.emplace_back(std::move(containedFace));
            status.increaseCounter(idSet3);
            status.print();
        }

        if (tryCount >= maxTriesForEntity) printSkipMessage();
    }

    // print the final entity density
    const auto domainVolume = computeMagnitude(domain);
    const auto numEntities = entitiesSet1.size()+entitiesSet2.size()+entitiesSet3.size();
    std::cout << "\nCreated a network consisting of " << numEntities << " entities\n";
    std::cout << "Volume of the domain: " << domainVolume << std::endl;
    std::cout << "Area of the network: " << area << " m²"
              << ", which corresponds to a density of " << area/domainVolume << " m²/m³" << std::endl;


    /////////////////////////////////////////////////////////////////////////////////////////
    // 5. The entities of the network have been created. We can now construct the network. //
    /////////////////////////////////////////////////////////////////////////////////////////

    // construct and write a contained network, i.e. write out both network and domain.
    std::cout << "Building and writing network" << std::endl;
    ContainedEntityNetworkBuilder<ctype> builder;

    // add sub-domains
    builder.addConfiningSubDomain(domain,     Id(1));
    builder.addSubDomainEntities(entitiesSet1, Id(1));
    builder.addSubDomainEntities(entitiesSet2, Id(1));
    builder.addSubDomainEntities(entitiesSet3, Id(1));

    // now we can build and write out the network in Gmsh file format
    GmshWriter gmshWriter(builder.build());
    gmshWriter.write("network", /*sizeAtEntities*/0.5, /*sizeInDomain*/5.0);

    // print time it took to generate the network
    const auto end = std::chrono::steady_clock::now();
    const auto duration = std::chrono::duration<double>(end-start).count();
    std::cout << "Overall CPU time was " << duration << " seconds." << std::endl;

    return 0;
}
