import frackit.geometry as geometry

# we want to use a cylindrical domain with radius=0.5 and height=1.0
domain = geometry.Cylinder(0.5, 1.0)

# we sample points uniformly within the shifted unit cube such
# that the origin is in the center of the bottom boundary of the domain
from frackit.sampling import makeUniformPointSampler
box = geometry.Box(-0.5, -0.5, 0.0, 0.5, 0.5, 1.0)
pointSampler = makeUniformPointSampler(box)

# returns a sampler from a gaussian distribution with mean and std deviation
def gaussianSampler(mean, stdDev):
    import random
    def sample():
        return random.gauss(mean, stdDev)
    return sample

#returns a sampler from a uniform distribution between min and max
def uniformSampler(min, max):
    import random
    def sample():
        return random.uniform(min, max)
    return sample

# we sample quadrialeterals within the box with gaussian distributions for the parameters
from frackit.common import toRadians
from frackit.sampling import QuadrilateralSampler as QuadSampler
quadSampler1 = QuadSampler(pointSampler,
                           gaussianSampler(toRadians(45.0), toRadians(5.0)), # strike angle
                           gaussianSampler(toRadians(90.0), toRadians(5.0)), # dip angle
                           uniformSampler(0.4, 0.6),                         # strike length
                           uniformSampler(0.4, 0.6))                         # dip length

# sampler for quadrilaterals of the secondary orientation
quadSampler2 = QuadSampler(pointSampler,
                           gaussianSampler(toRadians(0.0), toRadians(5.0)), # strike angle
                           gaussianSampler(toRadians(0.0), toRadians(5.0)), # dip angle
                           uniformSampler(0.4, 0.6),                        # strike length
                           uniformSampler(0.4, 0.6))                        # dip length

# constructs a constraints object with default settings for this example
from frackit.entitynetwork import EntityNetworkConstraints
def makeDefaultConstraints():
    c = EntityNetworkConstraints()
    c.setMinDistance(0.05)
    c.setMinIntersectingAngle(toRadians(30.0))
    c.setMinIntersectionMagnitude(0.05)
    c.setMinIntersectionDistance(0.05)
    return c

# We want to enforce some constraints on the set of quadrilaterals.
# In particular, for entities of the same set we want a minimum spacing
# distance of 5cm, and the quadrilaterals must not intersect in angles
# less than 30°. Moreover, if they intersect, we don't want intersection
# edges whose length is smaller than 5cm, and, the intersection should not
# be too close to the boundary of one of two intersecting quadrilaterals. Here: 5cm.
constraintsOnSelf = makeDefaultConstraints()

# with respect to entities of the other set, we want to have larger intersection angles
constraintsOnOther = makeDefaultConstraints()
constraintsOnOther.setMinIntersectingAngle(toRadians(40.0));

# we use the default constraints w.r.t. to the domain boundary
constraintsOnBoundary = makeDefaultConstraints()

# we use unique identifiers for the entities of the two orientations
from frackit.common import Id
idSet1 = Id(1)
idSet2 = Id(2)

# use the status class to define when to stop sampling
from frackit.sampling import SamplingStatus
status = SamplingStatus()
status.setTargetCount(idSet1, 15) # we want 15 entities in each set
status.setTargetCount(idSet2, 15) # we want 15 entities in each set

# start sampling into set 1 and keep alternating
sampleIntoSet1 = True

# lists to store the sampled entities
entities1 = []
entities2 = []

print("\n --- Start entity sampling ---\n")
while not status.finished():
    # sample a quadrilateral, alternating between sampler 1 and sampler 2
    quad = quadSampler1() if sampleIntoSet1 else quadSampler2()
    entitySet = entities1 if sampleIntoSet1 else entities2
    otherEntitySet = entities2 if sampleIntoSet1 else entities1

    # skip the rest if constraints are violated
    if not constraintsOnSelf.evaluate(entitySet, quad):
        status.increaseRejectedCounter()
        continue

    # skip the rest if constraints are violated
    if not constraintsOnOther.evaluate(otherEntitySet, quad):
        status.increaseRejectedCounter()
        continue

    # enforce constraints w.r.t. to the domain boundary
    if not constraintsOnBoundary.evaluate(domain.topFace(), quad):
        status.increaseRejectedCounter()
        continue
    if not constraintsOnBoundary.evaluate(domain.bottomFace(), quad):
        status.increaseRejectedCounter()
        continue

    from frackit.occutilities import getShape
    if not constraintsOnBoundary.evaluate(getShape(domain.lateralFace()), quad):
        status.increaseRejectedCounter()
        continue

    # reject entities whose contained part is too small
    from frackit.geometry import computeContainedMagnitude
    containedArea = computeContainedMagnitude(quad, domain);
    if (containedArea < 0.01):
        status.increaseRejectedCounter()
        continue

    # the entity is admissible
    entitySet.append(quad)
    id = idSet1 if sampleIntoSet1 else idSet2
    status.increaseCounter(id)
    status.print()

    # sample into other set the next time
    sampleIntoSet1 = not sampleIntoSet1

print("\n --- Finished entity sampling ---\n")

# We can now create an entity network from the two sets
from frackit.entitynetwork import ContainedEntityNetworkBuilder
builder = ContainedEntityNetworkBuilder()

# our domain confines all entities and receives the id 1
builder.addConfiningSubDomain(domain, Id(1));

# define all entities to be embedded in this domain
builder.addSubDomainEntities(entities1, Id(1));
builder.addSubDomainEntities(entities2, Id(1));

# let the builder construct the network and write it to gmsh file format
print("\n --- Constructing entity network from the raw entities ---\n")
network = builder.build();

print("\n --- Writing .geo file ---\n")
from frackit.io import GmshWriter
writer = GmshWriter(network);
writer.write("network", # filename of the .geo files (will add extension .geo automatically)
             0.1);      # element size to be used
print("\n --- Finished writing .geo file ---\n")
