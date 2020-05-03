from ._sampling import *


class BoxPointSampler:
    """Class to sample random points within a box."""

    def __init__(self, samplerX, samplerY, samplerZ):
        """Create the sampler from random variable samplers for the coordinate directions."""
        self.samplerX = samplerX
        self.samplerY = samplerY
        self.samplerZ = samplerZ

    def sample(self):
        x = self.samplerX()
        y = self.samplerY()
        z = self.samplerZ()

        from frackit.geometry import Point_3
        return Point_3(x, y, z)

# creator for a uniform sampler within an interval
def uniformIntervalSamplerCreator():
    def makeSampler(min, max):
        import random
        def doSample():
            return random.uniform(min, max)
        return doSample
    return makeSampler

def makeBoxPointSampler(box,
                        samplerCreatorX=uniformIntervalSamplerCreator(),
                        samplerCreatorY=uniformIntervalSamplerCreator(),
                        samplerCreatorZ=uniformIntervalSamplerCreator()):
    """Creates a BoxPointSampler using the provided creators for samplers on intervals.
       The creators default to uniform interval sampler creators if nothing is specified."""

    # make samplers
    samplerX = samplerCreatorX(box.xMin(), box.xMax())
    samplerY = samplerCreatorY(box.yMin(), box.yMax())
    samplerZ = samplerCreatorZ(box.zMin(), box.zMax())

    return BoxPointSampler(samplerX, samplerY, samplerZ)
