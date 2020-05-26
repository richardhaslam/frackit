import frackit.geometry as geometry

# we use the unit cube as domain
box = geometry.Box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)

# we sample points uniformly within the domain
from frackit.sampling import makeUniformPointSampler
pointSampler = makeUniformPointSampler(box)

# returns a sampler from a gaussian distribution with mean and std deviation
def gaussianSampler(mean, stdDev):
    import random
    def sample():
        return random.gauss(mean, stdDev)
    return sample

# we sample quadrialeterals within the box with gaussian distributions for the parameters
from frackit.common import toRadians
from frackit.sampling import QuadrilateralSampler as QuadSampler
quadSampler1 = QuadSampler(pointSampler,
                           gaussianSampler(toRadians(45.0), toRadians(5.0)), # strike angle
                           gaussianSampler(toRadians(90.0), toRadians(5.0)), # dip angle
                           gaussianSampler(0.5, 0.1),                        # edge length
                           0.05)                                             # threshold for minimum edge length

# sampler for quadrilaterals of the secondary orientation
quadSampler2 = QuadSampler(pointSampler,
                           gaussianSampler(toRadians(0.0), toRadians(5.0)), # strike angle
                           gaussianSampler(toRadians(0.0), toRadians(5.0)), # dip angle
                           gaussianSampler(0.5, 0.1),                       # edge length
                           0.05)                                            # threshold for minimum edge length

# We want to enforce some constraints on the set of quadrilaterals.
# In particular, for entities of the same set we want a minimum spacing
# distance of 5cm, and the quadrilaterals must not intersect in angles
# less than 30°. Moreover, if they intersect, we don't want intersection
# edges whose length is smaller than 5cm, and, the intersection should not
# be too close to the boundary of one of two intersecting quadrilaterals. Here: 5cm.
from frackit.entitynetwork import EntityNetworkConstraints
constraintsOnSelf = EntityNetworkConstraints()
constraintsOnSelf.setMinDistance(0.05)
constraintsOnSelf.setMinIntersectingAngle(toRadians(30.0))
constraintsOnSelf.setMinIntersectionMagnitude(0.05)
constraintsOnSelf.setMinIntersectionDistance(0.05)

# with respect to entities of the other set, we want to have larger intersection angles
constraintsOnOther = EntityNetworkConstraints()
constraintsOnOther.setMinDistance(0.05)
constraintsOnOther.setMinIntersectingAngle(toRadians(40.0));
constraintsOnOther.setMinIntersectionMagnitude(0.05)
constraintsOnOther.setMinIntersectionDistance(0.05)

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
    quad = quadSampler1.sample() if sampleIntoSet1 else quadSampler2.sample()
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

    # the entity is admissible
    entitySet.append(quad)
    id = idSet1 if sampleIntoSet1 else idSet2
    status.increaseCounter(id)
    status.print()

    # sample into other set the next time
    sampleIntoSet1 = not sampleIntoSet1

print("\n --- Finished entity sampling ---\n")

# We can now create an entity network from the two sets
from frackit.entitynetwork import EntityNetworkBuilder
builder = EntityNetworkBuilder()
builder.addEntities(entities1);
builder.addEntities(entities2);

# let the builder construct the network and write it to gmsh file format
print("\n --- Constructing entity network from the raw entities ---\n")
network = builder.build();

print("\n --- Writing .geo file ---\n")
from frackit.io import GmshWriter
writer = GmshWriter(network);
writer.write("network", # filename of the .geo files (will add extension .geo automatically)
             0.1);      # element size to be used
print("\n --- Finished writing .geo file ---\n")
