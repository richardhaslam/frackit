import sys
from time import process_time
startTime = process_time()

# print welcome message
print("\n\n"
      "#########################################################################\n"
      "## Example 4: Creation of polygonal entities around a tunnel structure ##\n"
      "#########################################################################\n"
      "\n\n")

###########################
# 1.Construct the domain ##
###########################

# In this example, we consider a fraction of a geological layer
# around a tunnel structure. To this end, we construct the geological
# layer as box, and intersect it with a sphere in order to restrict
# the model domain to a smaller region.
from frackit.geometry import Box, Sphere, Cylinder, Direction, Vector, Circle, Point

# the layer has a thickness of 40m, and initially we consider 100m lateral extent
layer = Box(-50.0, -50.0, -20.0, 50.0, 50.0, 20.0)
domainSphere = Sphere(Point(0.0, 0.0, 0.0), 25.0)

# a tunnel with a radius of 5m crosses the domain
tunnelRadius = 5.0
tunnelDirection = Direction(Vector(1.0, 0.0, 0.0))
tunnelBase = Circle(Point(-25.0, 0.0, 0.0), tunnelDirection, tunnelRadius)
tunnel = Cylinder(tunnelBase, 50.0)

# the domain is the layer constrained to the domainSphere and cut by the tunnel
from frackit.occutilities import getShape, getSolids, getFaces, cut
from frackit.occutilities import intersect as intersectShapes
solids = getSolids( intersectShapes(getShape(layer), getShape(domainSphere), 1e-6) )
if len(solids) != 1: sys.exit("Intersection operation with domainSphere should yield a single solid")

solids = getSolids( cut(solids[0], getShape(tunnel), 1e-6) )
if len(solids) != 1: sys.exit("Cut operation with tunnel should yield a single solid")

domain = solids[0]
tunnelSurfaceShape = getShape(tunnel.lateralFace())
sphereSurfaceShape = getFaces(getShape(domainSphere))
if len(sphereSurfaceShape) != 1: sys.exit("Sphere is expected to contain a single face")
sphereSurfaceShape = sphereSurfaceShape[0]

########################################################################
## 2. Make sampler classes to randomly sample points (entity centers) ##
##    and polygons as entities (using the sampled center points)      ##
########################################################################

# creates a uniform sampler within an interval
import random
def uniformSampler(a, b):
    def sample(): return random.uniform(a, b)
    return sample

# creates a sampler from a normal distribution with the given mean and std deviation
def normalSampler(mean, stdDev):
    def sample(): return random.gauss(mean, stdDev)
    return sample

# creates a sampler from an integer uniform distribution
def choiceSampler(a, b):
    def sample(): return random.choice([v for v in range(a, b+1)])
    return sample

# function to create polygon samplers for the different orientations
from frackit.common import toRadians
from frackit.sampling import PolygonSampler

sizeDistro = uniformSampler(3.0, 6.0)
numCornersDistro = choiceSampler(4, 9)

def getSampler(pointSampler, orientation):
    if orientation == 1:
        return PolygonSampler(pointSampler,                                    # sampler for polygon centers
                              normalSampler(toRadians(15.0), toRadians(10.0)), # strike angle: mean value & standard deviation
                              normalSampler(toRadians(90.0), toRadians(10.0)), # dip angle: mean value & standard deviation
                              sizeDistro, sizeDistro, numCornersDistro)        # strike & dip length, number of corners
    elif orientation == 2:
        return PolygonSampler(pointSampler,                                   # sampler for polygon centers
                              normalSampler(toRadians(0.0), toRadians(10.0)), # strike angle: mean value & standard deviation
                              normalSampler(toRadians(0.0), toRadians(10.0)), # dip angle: mean value & standard deviation
                              sizeDistro, sizeDistro, numCornersDistro)       # strike & dip length, number of corners
    elif orientation == 3:
        return PolygonSampler(pointSampler,                                     # sampler for polygon centers
                              normalSampler(toRadians(-75.0), toRadians(10.0)), # strike angle: mean value & standard deviation
                              normalSampler(toRadians(90.0), toRadians(10.0)),  # dip angle: mean value & standard deviation
                              sizeDistro, sizeDistro, numCornersDistro)         # strike & dip length, number of corners
    else:
        import sys
        sys.exit("Unsupported orientation index")


#######################################################################
## 3. Define constraints that should be fulfilled among the entities ##
##    of different orientations.                                     ##
#######################################################################

from frackit.entitynetwork import EntityNetworkConstraints, EntityNetworkConstraintsMatrix

# constructs default constraints for this test to be modified after
def makeDefaultConstraints():
    c = EntityNetworkConstraints()
    c.setMinDistance(1.0)
    c.setMinIntersectingAngle(toRadians(25.0))
    c.setMinIntersectionMagnitude(0.1)
    c.setMinIntersectionDistance(0.1)
    return c

# Define constraints between entities of orientation 1
constraints1 = makeDefaultConstraints()

# Define constraints between entities of orientation 2
# we want to enforce larger spacing between those entities
constraints2 = makeDefaultConstraints()
constraints2.setMinDistance(0.5)

# Define constraints between entities of orientation 3
constraints3 = makeDefaultConstraints()
constraints3.setMinDistance(0.15)

# Define constraints between entities of different sets
constraintsOnOther = makeDefaultConstraints()
constraintsOnOther.setMinDistance(0.25)
constraintsOnOther.setMinIntersectingAngle(toRadians(30.0))

# Moreover, we define constraints w.r.t. the domain boundary
constraintsOnBoundary = makeDefaultConstraints()
constraintsOnBoundary.setMinIntersectingAngle(toRadians(10.0))

# When polygons intersect, the intersection might be very close to
# one of the corner points. We want to avoid small length scales caused
# by this, and this lambda provides the check for it.
def producesSmallLengthScale(geometry, isection):

    from frackit.geometry import computeDistance
    from frackit.occutilities import getVertices
    verts = getVertices(getShape(geometry))
    if any((computeDistance(v, isection) < 1e-3 for v in verts)): return True
    return False

###########################
## 4. Network generation ##
###########################

# dictionary to store entities of the three orientations in
from frackit.common import Id
setId1 = Id(1); setId2 = Id(2); setId3 = Id(3);
entitySets = {setId1: [], setId2: [], setId3: []}

# In a first step, we want to populate entities of orientation 1 around
# the tunnel, and require each entity to be intersecting with it. To this end,
# we sample center points within a hollow cylinder around the tunnel.
from frackit.sampling import SamplingStatus
from frackit.geometry import HollowCylinder
numTargetEntities = 75
status = SamplingStatus()
status.setTargetCount(setId1, numTargetEntities)
sampleCylinder1 = HollowCylinder(tunnel.bottomFace().center(), tunnel.direction(),
                                 tunnel.radius(), tunnel.radius() + 2.0, tunnel.height())

from frackit.geometry import computeMagnitude, getBoundingBox
from frackit.geometry import intersect
from frackit.occutilities import getFaces
from frackit.sampling import makeUniformPointSampler
print('## Step 1: generation of entities of orientation 1 around the tunnel\n\n\n')

area = 0.0;
polygonSampler1 = getSampler(makeUniformPointSampler(sampleCylinder1), 1)
while not status.finished():
    candidate = polygonSampler1()
    isection = intersect(tunnel.lateralFace(), candidate)

    # skip entities that don't intersect the tunnel
    if all( "Empty" in g.name() for g in isection ):
        status.increaseRejectedCounter(); continue;

    # avoid small length scales produced by the intersections
    if producesSmallLengthScale(candidate, isection):
         status.increaseRejectedCounter(); continue;

    # compute the part contained in the domain
    containedShape = intersectShapes(getShape(candidate), domain, 1e-6)
    containedFaces = getFaces(containedShape)

    # skip everything that does not lead to a single contained face
    if len(containedFaces) != 1: status.increaseRejectedCounter(); continue;

    # if too little of the candidate is inside the domain, reject it
    containedFace = containedFaces[0];
    containedFaceArea = computeMagnitude(containedFace);
    if containedFaceArea < 0.25*candidate.area():
        status.increaseRejectedCounter(); continue;

    # check constraints to other entities and domain boundary
    if not constraints1.evaluate(entitySets[setId1], containedFace):
        status.increaseRejectedCounter(); continue;
    if not constraintsOnBoundary.evaluate(tunnelSurfaceShape, containedFace):
        status.increaseRejectedCounter(); continue;
    if not constraintsOnBoundary.evaluate(sphereSurfaceShape, containedFace):
        status.increaseRejectedCounter(); continue;

    area += containedFaceArea
    entitySets[setId1].append(containedFace)
    status.increaseCounter(setId1)
    status.print()

print('\n\n\n## Step 2: for each entity so far, try to generate an intersecting entity with orientation 2\n\n\n')

# after rejecting 500 candidates for an entity, we move to the next one
maxTriesForEntity = 500;
def printSkipMessage ():
    print('\t\t -- Skipped search for this entity after ' + \
          str(maxTriesForEntity) + ' unsuccessful tries. --\n')

from math import sqrt
status.reset()
status.setTargetCount(setId2, len(entitySets[setId1]))
for e in entitySets[setId1]:

    # sample an entity of orientation 2 in the vicinity of e
    rawBBox = getBoundingBox(e);
    charLength = sqrt(computeMagnitude(e));
    sampleBox = Box(rawBBox.xMin()-charLength, rawBBox.yMin()-charLength, rawBBox.zMin()-charLength,
                    rawBBox.xMax()+charLength, rawBBox.yMax()+charLength, rawBBox.zMax()+charLength);

    tryCount = 0;
    accepted = False;
    sampler = getSampler(makeUniformPointSampler(sampleBox), 2);

    while not accepted and tryCount <= maxTriesForEntity:
        tryCount += 1
        candidate = sampler()
        isection = intersect(candidate, e)

        # we only want intersecting entities
        if all( "Empty" in g.name() for g in isection ):
            status.increaseRejectedCounter(); continue;

        # avoid small length scales produced by the intersections
        if producesSmallLengthScale(candidate, isection):
             status.increaseRejectedCounter(); continue;
        if producesSmallLengthScale(e, isection):
             status.increaseRejectedCounter(); continue;

        # compute the part contained in the domain
        containedShape = intersectShapes(getShape(candidate), domain, 1e-6)
        containedFaces = getFaces(containedShape)

        # skip everything that does not lead to a single contained face
        if len(containedFaces) != 1: status.increaseRejectedCounter(); continue;

        # if too little of the candidate is inside the domain, reject it
        containedFace = containedFaces[0];
        containedFaceArea = computeMagnitude(containedFace);
        if containedFaceArea < 0.25*candidate.area():
            status.increaseRejectedCounter(); continue;

        # check constraints to other entities and domain boundary
        if not constraintsOnOther.evaluate(entitySets[setId1], containedFace):
            status.increaseRejectedCounter(); continue;
        if not constraints2.evaluate(entitySets[setId2], containedFace):
            status.increaseRejectedCounter(); continue;
        if not constraintsOnBoundary.evaluate(sphereSurfaceShape, containedFace):
            status.increaseRejectedCounter(); continue;
        if not constraintsOnBoundary.evaluate(tunnelSurfaceShape, containedFace):
            status.increaseRejectedCounter(); continue;

        accepted = True;
        area += containedFaceArea
        entitySets[setId2].append(containedFace)
        status.increaseCounter(setId2)
        status.print()

    if tryCount >= maxTriesForEntity: printSkipMessage()


print('\n\n\n## Step 3: for each entity of orientation 2, try to generate an intersecting entity with orientation 3\n\n\n')
status.reset()
status.setTargetCount(setId3, len(entitySets[setId2]))

for e  in entitySets[setId2]:

    # sample an entity of orientation 2 in the vicinity of e
    rawBBox = getBoundingBox(e);
    charLength = sqrt(computeMagnitude(e));
    sampleBox = Box(rawBBox.xMin()-charLength, rawBBox.yMin()-charLength, rawBBox.zMin()-charLength,
                    rawBBox.xMax()+charLength, rawBBox.yMax()+charLength, rawBBox.zMax()+charLength);

    tryCount = 0;
    accepted = False;
    sampler = getSampler(makeUniformPointSampler(sampleBox), 3);

    while not accepted and tryCount <= maxTriesForEntity:
        tryCount += 1
        candidate = sampler()
        isection = intersect(candidate, e)

        # we only want intersecting entities
        if all( "Empty" in g.name() for g in isection ):
            status.increaseRejectedCounter(); continue;

        # avoid small length scales produced by the intersections
        if producesSmallLengthScale(candidate, isection):
             status.increaseRejectedCounter(); continue;
        if producesSmallLengthScale(e, isection):
             status.increaseRejectedCounter(); continue;

        # compute the part contained in the domain
        containedShape = intersectShapes(getShape(candidate), domain, 1e-6)
        containedFaces = getFaces(containedShape)

        # skip everything that does not lead to a single contained face
        if len(containedFaces) != 1: status.increaseRejectedCounter(); continue;

        # if too little of the candidate is inside the domain, reject it
        containedFace = containedFaces[0];
        containedFaceArea = computeMagnitude(containedFace);
        if containedFaceArea < 0.25*candidate.area():
            status.increaseRejectedCounter(); continue;

        # check constraints to other entities and domain boundary
        if not constraintsOnOther.evaluate(entitySets[setId1], containedFace):
            status.increaseRejectedCounter(); continue;
        if not constraintsOnOther.evaluate(entitySets[setId2], containedFace):
            status.increaseRejectedCounter(); continue;
        if not constraints3.evaluate(entitySets[setId3], containedFace):
            status.increaseRejectedCounter(); continue;
        if not constraintsOnBoundary.evaluate(sphereSurfaceShape, containedFace):
            status.increaseRejectedCounter(); continue;
        if not constraintsOnBoundary.evaluate(tunnelSurfaceShape, containedFace):
            status.increaseRejectedCounter(); continue;

        accepted = True;
        area += containedFaceArea
        entitySets[setId3].append(containedFace)
        status.increaseCounter(setId3)
        status.print()

    if tryCount >= maxTriesForEntity: printSkipMessage()

# print the final entity density
domainVolume = computeMagnitude(domain);
numEntities = len(entitySets[setId1]) + len(entitySets[setId2]) + len(entitySets[setId3])
print("\nCreated a network consisting of {:d} entities.".format(numEntities))
print("Volume of the domain: {:.2f}.".format(domainVolume))
print("Area of the network: {:.2f} m², which corresponds to a density of {:.2f} m²/m³".format(area, area/domainVolume))

########################################################################################
## 5. The entities of the network have been created. We can now construct the network ##
########################################################################################

print("Building and writing network")
from frackit.entitynetwork import ContainedEntityNetworkBuilder
builder = ContainedEntityNetworkBuilder();

# add sub-domains
builder.addConfiningSubDomain(domain,    Id(1))
for id in entitySets: builder.addSubDomainEntities(entitySets[id], Id(1))

# now we can build and write out the network in Gmsh file format
from frackit.io import GmshWriter
gmshWriter = GmshWriter(builder.build());
gmshWriter.write("network", 0.5, 5.0)

stopTime = process_time()
print("Overall CPU time was {:.2f} seconds.".format(stopTime-startTime))
