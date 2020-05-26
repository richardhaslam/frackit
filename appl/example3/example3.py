# print welcome message
print("\n\n"
      "##############################################################################\n"
      "## Example: Creation of entities (disks/quadrilaterals) in a layered medium ##\n"
      "##############################################################################\n"
      "\n\n")

####################################################
# 1. Read in the domain geometry from .brep file. ##
#    The file name is defined in CMakeLists.txt   ##
####################################################
from frackit.occutilities import readShape, getSolids, getShells, getFaces, getBoundingBox
domainShape = readShape("layers.brep")

# obtain the three solids contained in the file
solids = getSolids(domainShape)
if len(solids) != 3: raise RuntimeError("Expected the .brep file to contain 3 solids")

# The sub-domain we want to create a network in is the center one.
# Compute its volume and get the boundary faces for constraint evaluation.
from frackit.magnitude import computeMagnitude
networkDomain = solids[1]
domainVolume = computeMagnitude(networkDomain);
shells = getShells(networkDomain)
if len(shells) != 1: raise RuntimeError("Expected a single shell bounding the domain")

domainBoundaryFaces = getFaces(shells[0])
if len(domainBoundaryFaces) != 6: raise RuntimeError("Expected 6 faces to bound the domain")


##############################################################################
## 2. Make sampler classes to randomly sample points (entity centers)       ##
##    and disk/quadrilaterals as entities (using the sampled center points) ##
##############################################################################

# Bounding box of the domain in which we want to place the entities
# In order for this to work we parse the solid back into a wrapper around a BRep shape
from frackit.geometry import Box
from frackit.occutilities import OCCShapeWrapper
domainBBox = getBoundingBox( OCCShapeWrapper(networkDomain) )

# creates a uniform sampler within an interval
import random
def uniformSampler(a, b):
    def sample(): return random.uniform(a, b)
    return sample

# creates a sampler from a normal distribution with the given mean and std deviation
def normalSampler(mean, stdDev):
    def sample(): return random.gauss(mean, stdDev)
    return sample

from frackit.common import toRadians, Id
from frackit.sampling import DiskSampler, QuadrilateralSampler, makeUniformPointSampler
# sampler for disks (orientation 1)
diskSampler = DiskSampler(makeUniformPointSampler(domainBBox),           # sampler for disk center points
                          uniformSampler(30.0, 6.5),                     # major axis length: mean value & standard deviation
                          uniformSampler(24.0, 4.5),                     # minor axis length: mean value & standard deviation
                          normalSampler(toRadians(0.0), toRadians(7.5)), # rotation around x-axis: mean value & standard deviation
                          normalSampler(toRadians(0.0), toRadians(7.5)), # rotation around y-axis: mean value & standard deviation
                          normalSampler(toRadians(0.0), toRadians(7.5))) # rotation around z-axis: mean value & standard deviation

# sampler for quadrilaterals (orientation 2)
quadSampler = QuadrilateralSampler(makeUniformPointSampler(domainBBox),             # sampler for quadrilateral center points
                                   normalSampler(toRadians(45.0), toRadians(5.0)),  # strike angle: mean value & standard deviation
                                   normalSampler(toRadians(90.0), toRadians(5.0)),  # dip angle: mean value & standard deviation
                                   normalSampler(45.0, 5.0),                        # edge length: mean value & standard deviation
                                   5.0)                                             # threshold for minimum edge length

# Define ids for the two entity sets
diskSetId = Id(1) # we give the set of orientation one, consisting of disks, the id 1
quadSetId = Id(2) # we give the set of orientation two, consisting of quadrilaterals, the id 2


#######################################################################
## 3. Define constraints that should be fulfilled among the entities ##
##    of different orientations.                                     ##
#######################################################################

from frackit.entitynetwork import EntityNetworkConstraints, EntityNetworkConstraintsMatrix

# constructs default constraints for this test to be modified after
def makeDefaultConstraints():
    c = EntityNetworkConstraints()
    c.setMinDistance(2.5)
    c.setMinIntersectingAngle(toRadians(25.0))
    c.setMinIntersectionMagnitude(5.0)
    c.setMinIntersectionDistance(2.5)
    return c

# Define constraints between entities of orientation 1
constraints1 = makeDefaultConstraints()

# Define constraints between entities of orientation 2
# we want to enforce larger spacing between those entities
constraints2 = makeDefaultConstraints()
constraints2.setMinDistance(5.0)

# Define constraints between entities of different sets
constraintsOnOther = makeDefaultConstraints()
constraintsOnOther.setMinDistance(2.5)
constraintsOnOther.setMinIntersectingAngle(toRadians(40.0))

# We can use the constraints matrix to facilitate constraint evaluation
constraintsMatrix = EntityNetworkConstraintsMatrix()
constraintsMatrix.addConstraints(constraints1,                       # constraint instance
                                 [diskSetId.get(), diskSetId.get()]) # sets between which to use these constraints

constraintsMatrix.addConstraints(constraints2,                       # constraint instance
                                 [quadSetId.get(), quadSetId.get()]) # sets between which to use these constraints

constraintsMatrix.addConstraints(constraintsOnOther,                   # constraint instance
                                 [[diskSetId.get(), quadSetId.get()],  # sets between which to use these constraints
                                  [quadSetId.get(), diskSetId.get()]]) # sets between which to use these constraints

# Moreover, we define constraints w.r.t. the domain boundary
constraintsOnDomain = makeDefaultConstraints()
constraintsOnDomain.setMinIntersectingAngle(toRadians(15.0))


###########################
## 4. Network generation ##
###########################

# Helper class for terminal output of the creation
# progress and definition of stop criterion etc
from frackit.sampling import SamplingStatus
status = SamplingStatus()
status.setTargetCount(diskSetId, 12) # we want 11 entities of orientation 1
status.setTargetCount(quadSetId, 16) # we want 13 entities of orientation 2

# store all entity sets in a dictionary
entitySets = {diskSetId.get(): [], quadSetId.get(): []}

# Alternate between set 1 & set 2 during sampling phase
sampleIntoSet1 = True
containedNetworkArea = 0.0;
while not status.finished():
    id = diskSetId if sampleIntoSet1 else quadSetId
    geom = diskSampler.sample() if sampleIntoSet1 else quadSampler.sample()

    # If the set this geometry belongs to is finished, skip the rest
    if status.finished(id):
        sampleIntoSet1 = not sampleIntoSet1
        status.increaseRejectedCounter()
        continue

    # Moreover, we want to avoid small fragments (< 250 m²)
    from frackit.magnitude import computeContainedMagnitude
    containedArea = computeContainedMagnitude(geom, networkDomain);
    if containedArea < 350.0:
        status.increaseRejectedCounter()
        continue

    # enforce constraints w.r.t. to the other entities
    if not constraintsMatrix.evaluate(entitySets, geom, id):
        status.increaseRejectedCounter()
        continue

    # enforce constraints w.r.t. the domain boundaries
    if not constraintsOnDomain.evaluate(domainBoundaryFaces, geom):
        status.increaseRejectedCounter()
        continue

    # the geometry is admissible
    entitySets[id.get()].append(geom)
    status.increaseCounter(id)
    status.print()

    # keep track of entity area
    containedNetworkArea += containedArea;

    # sample from the other set next time
    sampleIntoSet1 = not sampleIntoSet1

# print the final entity density
density = containedNetworkArea/domainVolume;
print("\nEntity density of the contained network: {:f} m²/m³\n".format(density))


##########################################################################
## 5. The entities of the network have been created. We can now         ##
##    construct different types of networks (contained, confined, etc.) ##
##########################################################################

# construct and write a contained network, i.e. write out both network and domain.
print("Building and writing contained, confined network\n")
from frackit.entitynetwork import EntityNetworkBuilder, ContainedEntityNetworkBuilder
builder = ContainedEntityNetworkBuilder();

# add sub-domains
builder.addConfiningSubDomain(solids[0],     Id(1));
builder.addConfiningSubDomain(networkDomain, Id(2));
builder.addConfiningSubDomain(solids[2],     Id(3));

# The entites that we created all belong to subdomain 2
for setId in entitySets: builder.addSubDomainEntities(entitySets[setId], Id(2))

# now we can build and write out the network in Gmsh file format
from frackit.io import GmshWriter
gmshWriter = GmshWriter(builder.build());
gmshWriter.write("contained_confined", # body of the filename to be used (will add .geo)
                 2.5,                  # mesh size to be used on entities
                 5.0);                 # mesh size to be used on domain boundaries

# we can also not confine the network to its sub-domain,
# simply by adding the sub-domains as non-confining
print("Building and writing contained, unconfined network\n")
builder.clear();
builder.addSubDomain(solids[0],     Id(1));
builder.addSubDomain(networkDomain, Id(2));
builder.addSubDomain(solids[2],     Id(3));
for setId in entitySets: builder.addSubDomainEntities(entitySets[setId], Id(2))

gmshWriter = GmshWriter(builder.build());
gmshWriter.write("contained_unconfined", 2.5, 5.0);

# We could also only write out the network, without the domain
# For example, confining the network to the sub-domain...
print("Building and writing uncontained, confined network")
uncontainedBuilder = EntityNetworkBuilder()
uncontainedBuilder.addConfiningSubDomain(networkDomain, Id(2));
for setId in entitySets: uncontainedBuilder.addSubDomainEntities(entitySets[setId], Id(2))

gmshWriter = GmshWriter(uncontainedBuilder.build());
gmshWriter.write("uncontained_confined", 2.5);

# ... or not confining it
print("Building and writing uncontained, unconfined network\n")
uncontainedBuilder.clear();
for setId in entitySets: uncontainedBuilder.addSubDomainEntities(entitySets[setId], Id(2))

gmshWriter = GmshWriter(uncontainedBuilder.build());
gmshWriter.write("uncontained_unconfined", 2.5);
