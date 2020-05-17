from ._sampling import *


class BoxPointSampler:
    """Class to sample random points within a box."""

    def __init__(self, samplerX, samplerY, samplerZ):
        """
        Create the sampler from random variable samplers for the coordinate directions.

        Parameters:
        samplerX: sampler that samples from a distribution for the x-coordinate
        samplerY: sampler that samples from a distribution for the y-coordinate
        samplerZ: sampler that samples from a distribution for the z-coordinate
        """
        self.samplerX = samplerX
        self.samplerY = samplerY
        self.samplerZ = samplerZ

    def sample(self):
        x = self.samplerX()
        y = self.samplerY()
        z = self.samplerZ()

        from frackit.geometry import Point_3
        return Point_3(x, y, z)


class CylinderPointSampler:
    """Class to sample random points within a cylinder."""

    def __init__(self, cylinder, samplerR2, samplerPhi, samplerZ):
        """
        Create the sampler from random variable samplers for the coordinate directions.

        Points on the cylinder are represented in cylinder coordinates, i.e. by a radial,
        an angular and a height coordinate. Points on the cylinder are sampled based on
        the given distributions for the squared radius, the angle and the height.

        Parameters:
        cylinder: an instance of a cylinder class
        samplerR2: sampler that samples from a distribution for the squared radius
        samplerPhi: sampler that samples from a distribution for the angular coordinate
        samplerZ: sampler that samples from a distribution for the height
        """
        self.cylinder = cylinder
        self.samplerR2 = samplerR2
        self.samplerPhi = samplerPhi
        self.samplerZ = samplerZ

    def sample(self):

        from frackit.geometry import Vector_3
        a = Vector_3(self.cylinder.bottomFace().majorAxis());
        b = Vector_3(self.cylinder.bottomFace().minorAxis());
        n = Vector_3(self.cylinder.bottomFace().normal());

        from math import sqrt, sin, cos
        r = sqrt(self.samplerR2())
        phi = self.samplerPhi()
        z = self.samplerZ()

        a *= r*cos(phi);
        b *= r*sin(phi);
        n *= z;

        result = self.cylinder.bottomFace().center();
        result += a;
        result += b;
        result += n;

        return result;


# creator for a uniform sampler within an interval
def uniformIntervalSamplerCreator():
    def makeSampler(min, max):
        import random
        def doSample():
            return random.uniform(min, max)
        return doSample
    return makeSampler


def makePointSampler(geometry,
                     samplerCreatorBase1=uniformIntervalSamplerCreator(),
                     samplerCreatorBase2=uniformIntervalSamplerCreator(),
                     samplerCreatorBase3=uniformIntervalSamplerCreator()):
    """
    Creates a point sampler for the provided three-dimensional geometry using the
    provided creators for samplers on intervals.

    Points on the three-dimensional geometries considered here can be expressed in
    terms of tuples of coordinates (a1, a2, a3), where the entries of the tuple refer
    to the coordinate values with respect to the geometry-local basis. Moreover,
    the considered geometries are finite, which is expressed in terms of lower and
    upper bounds for the coordinates. For instance, the points within an axis-aligned
    box fall into the intervals ([xMin, xMax], [yMin, yMax], [zMin, zMax]). For a
    cylinder, described in cylinder coordinates, all points fall into the intervals
    ([0, r], [0, 2*Pi], [0, H]), where r is the radius an H is the height of the
    cylinder. The provided sampler creators are used to build samplers that sample
    coordinates within these intervals following a user-defined distribution. For
    more details on the meaning of the samplers, have a look at the geometry-specific
    point sampler implementations.

    The creators default to uniform interval sampler creators if nothing is specified.

    Parameters:
    geometry: the geometry for which a point sampler is to be created
    samplerCreatorBase1: creator for a sampler of coordinate values for the first coordinate
    samplerCreatorBase2: creator for a sampler of coordinates values for the second coordinate
    samplerCreatorBase3: creator for a sampler of coordinates values for the third coordinate
    """

    # box point sampler
    if geometry.name() == "Box":
        samplerX = samplerCreatorBase1(geometry.xMin(), geometry.xMax())
        samplerY = samplerCreatorBase2(geometry.yMin(), geometry.yMax())
        samplerZ = samplerCreatorBase3(geometry.zMin(), geometry.zMax())
        return BoxPointSampler(samplerX, samplerY, samplerZ)

    # cylinder point sampler
    if geometry.name() == "Cylinder":
        import math
        samplerR2  = samplerCreatorBase1(0.0, geometry.radius())
        samplerPhi = samplerCreatorBase2(0.0, 2.0*math.pi)
        samplerZ   = samplerCreatorBase3(0.0, geometry.height())
        return CylinderPointSampler(geometry, samplerR2, samplerPhi, samplerZ)

    # hollow cylinder point sampler
    if geometry.name() == "HollowCylinder":
        import math
        samplerR2  = samplerCreatorBase1(geometry.innerRadius(), geometry.outerRadius())
        samplerPhi = samplerCreatorBase2(0.0, 2.0*math.pi)
        samplerZ   = samplerCreatorBase3(0.0, geometry.height())
        return CylinderPointSampler(geometry.fullCylinder(), samplerR2, samplerPhi, samplerZ)

    raise NotImplementedError("No point sample creation formula implemented for provided geometry: " + geometry.name())
