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

    def __call__(self):
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

    def __call__(self):

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



class DiskSampler:
    """Class to randomly sample disk geometries."""

    def __init__(self, pointSampler, majAxisSampler, minAxisSampler, xAngleSampler, yAngleSampler, zAngleSampler):
        """
        Create the sampler from random variable samplers for the geometric properties.

        Parameters:
        pointSampler: sampler for points to be used as the disk centers
        majAxisSampler: samples from a distribution for the major axis length
        minAxisSampler: samples from distribution the minor axis length
        xAngleSampler: samples from a distribution for the angle or rotation around x-axis
        yAngleSampler: samples from a distribution for the angle or rotation around y-axis
        zAngleSampler: samples from a distribution for the angle or rotation around z-axis
        """
        self.pointSampler = pointSampler
        self.majAxisSampler = majAxisSampler
        self.minAxisSampler = minAxisSampler
        self.xAngleSampler = xAngleSampler
        self.yAngleSampler = yAngleSampler
        self.zAngleSampler = zAngleSampler

    def __call__(self):

        a = self.majAxisSampler()
        while (a <= 0.0): a = self.majAxisSampler()

        b = self.minAxisSampler()
        while (b <= 0.0): b = self.minAxisSampler()

        if (b > a): b = a

        alpha = self.xAngleSampler()
        beta = self.yAngleSampler()
        gamma = self.zAngleSampler()

        # find major/minor axis by rotations
        from frackit.geometry import Vector_3
        axes = [Vector_3(1.0, 0.0, 0.0), Vector_3(0.0, 1.0, 0.0)]

        from frackit.geometry import Direction_3
        e1 = Direction_3(axes[0])
        e2 = Direction_3(axes[1])
        e3 = Direction_3(Vector_3(0.0, 0.0, 1.0))

        from frackit.common import rotate
        rotate(axes[1], e1, alpha)  # rotate minor axis around x
        rotate(axes[0], e2, beta)   # rotate both axes around y
        rotate(axes[1], e2, beta)   # rotate both axes around y
        rotate(axes[0], e3, gamma)  # rotate both axes around z
        rotate(axes[1], e3, gamma)  # rotate both axes around z

        # sample center point and make disk
        from frackit.geometry import Ellipse_3
        ellipse = Ellipse_3(self.pointSampler(),
                            Direction_3(axes[0]),
                            Direction_3(axes[1]),
                            a, b)
        from frackit.geometry import Disk
        return Disk(ellipse)

class PolygonSampler:
    """Class to randomly sample polygons in 3d space."""

    def __init__(self, pointSampler, strikeAngleSampler, dipAngleSampler, strikeLengthSampler, dipLengthSampler, numCornersSampler):
        """
        Create the sampler from random variable samplers for the geometric properties.

        Parameters:
        pointSampler: sampler for points to be used as the qudrilateral centers
        strikeAngleSampler: samples from a distribution for the strike angle
        dipAngleSampler: samples from a distribution for the dip angle
        strikeLengthSampler: samples from distribution the dimensions in strike direction
        dipLengthSampler: samples from distribution the dimensions in strike direction
        numCornersSampler: samples from distribution the number of polygon corners
        """
        self.pointSampler = pointSampler
        self.strikeAngleSampler = strikeAngleSampler
        self.dipAngleSampler = dipAngleSampler
        self.strikeLengthSampler = strikeLengthSampler
        self.dipLengthSampler = dipLengthSampler
        self.numCornersSampler = numCornersSampler
        self.intervalSampler = uniformIntervalSamplerCreator()(0.0, 1.0)

    def __call__(self):
        from frackit.geometry import Polygon_3
        return Polygon_3(self.sampleCorners())

    def sampleCorners(self):
        count = 0
        numCorners = 0
        while numCorners < 3 and count < 100:
            numCorners = self.numCornersSampler()

        if count == 100:
            raise RuntimeError("Did not obtain more than 3 corners after sampling " \
                               "100 times from the provided distribution.")

        # randomly sample within interval [0, 1]
        xValues = [self.intervalSampler() for i in range(0, numCorners)]
        yValues = [self.intervalSampler() for i in range(0, numCorners)]

        # make sure interval is exploited
        xMin = min(xValues); xValues.remove(xMin); xMin = 0.0 if xMin > 0.1 else xMin
        yMin = min(yValues); yValues.remove(yMin); yMin = 0.0 if yMin > 0.1 else yMin
        xMax = max(xValues); xValues.remove(xMax); xMax = 1.0 if xMax < 0.9 else xMax
        yMax = max(yValues); yValues.remove(yMax); yMax = 1.0 if yMax < 0.9 else yMax

        # split values into two lists
        xVals = [[xMin], [xMin]]
        yVals = [[yMin], [yMin]]
        for i in range(0, numCorners-2):
            xVals[i%2].append(xValues[i])
            yVals[i%2].append(yValues[i])

        # sort the arrays and add max values
        for i in [0, 1]: xVals[i].sort(); xVals[i].append(xMax)
        for i in [0, 1]: yVals[i].sort(); yVals[i].append(yMax)

        # transform values into deltas
        for vals in [xVals[0], xVals[1], yVals[0], yVals[1]]:
            for i in range(0, len(vals)-1):
                vals[i] = vals[i+1] - vals[i]
            vals.pop()

        # flip sign of second arrays
        xVals[1] = [ -1.0*val for val in xVals[1] ]
        yVals[1] = [ -1.0*val for val in yVals[1] ]

        # make two arrays of delta values
        deltaX = [dx for dx in xVals[0]]; deltaX.extend(xVals[1])
        deltaY = [dy for dy in yVals[0]]; deltaY.extend(yVals[1])

        # compute the maximum deltaX/deltaY values
        deltaX.sort()
        maxDeltaX = max(deltaX)
        maxDeltaY = max(deltaY)
        maxDx2 = maxDeltaX*maxDeltaX + maxDeltaY*maxDeltaY

        deltas = []
        for dx in deltaX:
            tmp = maxDx2 - dx*dx

            # look for dy such that dx^2+dy^2 closest to maxDx2
            closest = maxDx2*100
            for dy in deltaY:
                diff = abs(tmp - dy*dy)
                if diff < closest:
                    closest = diff
                    closestDy = dy

            deltas.append([dx, closestDy])
            deltaY.remove(closestDy)

        if len(deltas) != numCorners:
            raise RuntimeError("Error during edge vectors construction")

        # scale vectors with actual dimensions (verify that they are > 0.0)
        strikeLength  = self.strikeLengthSampler()
        dipLength = self.dipLengthSampler()

        countStrike = 0;
        while strikeLength < 0.0 and countStrike < 100:
            strikeLength = self.strikeLengthSampler()
            countStrike += 1
        if countStrike == 100: raise RuntimeError("Could not sample strike dimension > 0.0 in 100 tries")

        countDip = 0;
        while strikeLength < 0.0 and countDip < 100:
            dipLength = self.dipLengthSampler()
            countDip += 1
        if countDip == 100: raise RuntimeError("Could not sample dip dimension > 0.0 in 100 tries")

        # construct edge vectors
        from frackit.geometry import Vector
        edgeVectors = [ Vector(dx*dipLength, dy*strikeLength, 0.0) for dx, dy in deltas ]

        # returns the angle w.r.t the standard x-y basis of an edge vector
        import math
        def getAngle(v):
            if abs(v.y()) < 1e-6: return 0.0 if v.x() > 0.0 else math.pi
            if abs(v.x()) < 1e-6: return math.pi/2.0 if v.y() > 0.0 else 3.0*math.pi/2.0

            angle = math.atan2(v.y(), v.x())
            return 2.0*math.pi + angle if v.y() < 0.0 else angle

        # sort edges by angle
        edgeVectors.sort(key=lambda v: getAngle(v))

        # rotate by strike and dip angle
        from frackit.common import rotate
        from frackit.geometry import Direction
        dipAngle = self.dipAngleSampler()
        strikeAngle = self.strikeAngleSampler()
        strikeDir = Direction(Vector(-math.sin(strikeAngle), math.cos(strikeAngle), 0.0))
        for v in edgeVectors: rotate(v, Direction(Vector(0.0, 0.0, 1.0)), strikeAngle)
        for v in edgeVectors: rotate(v, strikeDir, dipAngle)

        # combine vectors to form a polygon
        from frackit.geometry import Point
        vertices = [Point(0.0, 0.0, 0.0)]

        # first vertex inserted above is origin
        xMin = 0.0; xMax = 0.0; yMin = 0.0; yMax = 0.0; zMin = 0.0; zMax = 0.0;
        for i in range(0, len(edgeVectors)-1):
            vertices.append(vertices[-1] + edgeVectors[i])
            xMin = min(xMin, vertices[-1].x()); xMax = max(xMax, vertices[-1].x());
            yMin = min(yMin, vertices[-1].y()); yMax = max(yMax, vertices[-1].y());
            zMin = min(zMin, vertices[-1].z()); zMax = max(zMax, vertices[-1].z());

        # make sure next and last vertices are identical
        if not (vertices[-1] + edgeVectors[-1]).isEqual(vertices[0]):
            print("First vertex: ", vertices[0])
            print("Last vertex: ", vertices[-1])
            print("Past last vertex: ", vertices[-1] + edgeVectors[-1])
            raise RuntimeError("Past last and first vertex of polygon are not equal!")

        # move the polygon to sampled center
        c = self.pointSampler()
        bboxCenter = Point(xMin + 0.5*(xMax - xMin),
                           yMin + 0.5*(yMax - yMin),
                           zMin + 0.5*(zMax - zMin))
        dispVector = Vector(bboxCenter, c)
        for v in vertices: v += dispVector

        return vertices

class QuadrilateralSampler:
    """Class to randomly sample quadrilaterals in 3d space."""

    from frackit.precision import Precision
    def __init__(self, pointSampler, strikeAngleSampler, dipAngleSampler, strikeLengthSampler, dipLengthSampler = Precision.confusion):
        """
        Create the sampler from random variable samplers for the geometric properties.
        Internally, an instance of PolygonSampler is used, to which the arguments are forwarded.
        Note: This constructor is still compatible with the deprecated quadrilateral sampler
              algorithm and its constructor arguments. A warning will be printed if it is used as such.

        Parameters:
        pointSampler: sampler for points to be used as the qudrilateral centers
        strikeAngleSampler: samples from a distribution for the strike angle
        dipAngleSampler: samples from a distribution for the dip angle
        strikeLengthSampler: samples from distribution the dimensions in strike direction
        dipLengthSampler: samples from distribution the dimensions in strike direction
        """

        # if the last argument is not callable, we assume one is using the deprecated
        # way of construction and algorithm, where the last argument is the minumum edge length
        try: tmp = dipLengthSampler()
        except:
            import warnings
            warnings.warn("Please use new quadrilateral sampler implementation", DeprecationWarning, stacklevel=2)
            self.pointSampler = pointSampler
            self.strikeAngleSampler = strikeAngleSampler
            self.dipAngleSampler = dipAngleSampler
            self.edgeLengthSampler = strikeLengthSampler
            self.minEdgeLength = dipLengthSampler
            self.useDeprecated = True
        else:
            import random
            def cornerSampler(): return 4
            self.polygonSampler = PolygonSampler(pointSampler,
                                                 strikeAngleSampler, dipAngleSampler,
                                                 strikeLengthSampler, dipLengthSampler,
                                                 cornerSampler)
            self.useDeprecated = False

    def __call__(self):

        if not self.useDeprecated:
            c = self.polygonSampler.sampleCorners()
            if len(c) != 4:
                raise RuntimeError("Size mismatch of corner list!")

            from frackit.geometry import Quadrilateral_3
            return Quadrilateral_3(c[0], c[1], c[3], c[2])
        else:
            strike = self.strikeAngleSampler()
            dip = self.dipAngleSampler()

            # get the basis of the plane within the x-y-plane by
            # rotation around the z-axis with the strike angle and
            # rotation of the resulting first axis around the second by dip angle
            from frackit.geometry import Vector_3
            axes = [Vector_3(1.0, 0.0, 0.0), Vector_3(0.0, 1.0, 0.0)]

            from frackit.geometry import Direction_3
            from frackit.common import rotate
            rotate(axes[0], Direction_3(Vector_3(0.0, 0.0, 1.0)), strike);
            rotate(axes[1], Direction_3(Vector_3(0.0, 0.0, 1.0)), strike);
            rotate(axes[0], Direction_3(axes[1]), dip);

            # sample edge lengths until all are admissible
            dx1 = self.edgeLengthSampler()
            dx2 = self.edgeLengthSampler()
            dy1 = self.edgeLengthSampler()
            dy2 = self.edgeLengthSampler()
            while (dx1 < self.minEdgeLength): dx1 = self.edgeLengthSampler()
            while (dx2 < self.minEdgeLength): dx2 = self.edgeLengthSampler()
            while (dy1 < self.minEdgeLength): dy1 = self.edgeLengthSampler()
            while (dy2 < self.minEdgeLength): dy2 = self.edgeLengthSampler()

            from copy import deepcopy
            dxVec1 = Vector_3(deepcopy(axes[0].x()), deepcopy(axes[0].y()), deepcopy(axes[0].z())); dxVec1 *= dx1/2.0
            dxVec2 = Vector_3(deepcopy(axes[0].x()), deepcopy(axes[0].y()), deepcopy(axes[0].z())); dxVec2 *= dx2/2.0
            dyVec1 = Vector_3(deepcopy(axes[1].x()), deepcopy(axes[1].y()), deepcopy(axes[1].z())); dyVec1 *= dy1/2.0
            dyVec2 = Vector_3(deepcopy(axes[1].x()), deepcopy(axes[1].y()), deepcopy(axes[1].z())); dyVec2 *= dy2/2.0

            # compute corner points
            c = self.pointSampler()
            from frackit.geometry import Point_3
            c1 = Point_3(deepcopy(c.x()), deepcopy(c.y()), deepcopy(c.z())); c1 -= dxVec1; c1 -= dyVec1
            c2 = Point_3(deepcopy(c.x()), deepcopy(c.y()), deepcopy(c.z())); c2 += dxVec2; c2 -= dyVec2
            c3 = Point_3(deepcopy(c.x()), deepcopy(c.y()), deepcopy(c.z())); c3 -= dxVec1; c3 += dyVec1
            c4 = Point_3(deepcopy(c.x()), deepcopy(c.y()), deepcopy(c.z())); c4 += dxVec2; c4 += dyVec2

            from frackit.geometry import Quadrilateral_3
            return Quadrilateral_3(c1, c2, c3, c4);
