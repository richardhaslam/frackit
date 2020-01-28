<!--- Example picture --->
<p align="center">
    <img src="../../doc/img/example2_network.png" alt="frackit example 1" width="800"/>
</p>

Example 2
=========

In this exemplary application, a network of quadrilaterals, following the same
distributions as the network of [example 1][0], is created. However, in this
example we want to confined the network to a cylindrical domain (see image above).
We represent the domain via an instance of the internal cylinder class:

```cpp
// In this example we consider a cylindrical domain
Cylinder<ctype> domain(/*radius*/0.5, /*height*/1.0);
```

This creates an axis-aligned cylinder where the x- and y- axis form the basis of
the bottom face of the cylinder, while the cylinder axis is aligned with the
z-axis. Note that there are other ways to construct arbitrarily-oriented cylinders,
and for details on this we refer to the [class documentation][2].

In addition to constraints among entities of the different orientations, we also
want to enforce constraints with respect to the boundary of the cylindrical domain.
In this example, we simply reuse the constraints of entities of the same orientation
and write

```cpp
auto constraintsOnBoundary = constraintsOnSelf;
```

to create a new constraints object. Within the loop in which the entities are created
these constraints are enforced:

```cpp
// we want to enforce constraints also w.r.t. to the cylinder boundary
if (!constraintsOnBoundary.evaluate(domain.topFace(), quad))
{ status.increaseRejectedCounter(); continue; }
if (!constraintsOnBoundary.evaluate(domain.bottomFace(), quad))
{ status.increaseRejectedCounter(); continue; }
if (!constraintsOnBoundary.evaluate(domain.lateralFace(), quad))
{ status.increaseRejectedCounter(); continue; }
```

As you can see, the `Cylinder` class provides functions for obtaining the representations
of the top, bottom and lateral boundaries. The top and bottom boundaries are represented
by instances of the `Disk` class, while the lateral surface is described by the class
`CylinderSurface`. Moreover, we want to reject all those quadrilaterals of which only
a very small portion is inside the cylinder. This is done in the lines

```cpp
const auto containedArea = computeContainedMagnitude(quad, domain);

// reject if this is too small (< 0.01 m^2)
if (containedArea < 0.01)
{ status.increaseRejectedCounter(); continue; }
```

where `computeContainedMagnitude()` is a free function that returns the length/area/volume
of the part of a geometry that is contained inside of another geometry.

Finally, after the desired number of entities has been created, we cast the entities
into an instance of a `ContainedEntityNetwork`, using the corresponding builder class:

```cpp
ContainedEntityNetworkBuilder builder;

// define the domain (single sub-domain) and give it a unique id
builder.addConfiningSubDomain(domain, Id(1));
```

The class `ContainedEntityNetworkBuilder` returns an instance of a `ContainedEntityNetwork`
when the `build()` function is called (see below). This network implementation contains information
on (sub-)domains and which entities are embedded in which (sub-)domain. Each (sub-)domain
receives a unique identifier by means of the `Id` class. By calling
`builder.addConfiningSubDomain(domain, Id(1));`, we define the (sub-)domain to be
confining, i.e. entities that are added to this (sub-)domain will be confined to it
by cutting away all parts that lie outside the (sub-)domain. In contrast to that,
one could call `builder.addSubDomain(domain, Id(1));`, in which case embedded networks
will not be confined (see [example 3][2]). Entities are associated with the (sub-)
domain they are embedded in, and are added to the builder class by writing

```cpp
// define entities to be embedded in this domain
builder.addSubDomainEntities(entitySet1, Id(1));
builder.addSubDomainEntities(entitySet2, Id(1));

const auto network = builder.build();
```

The variable `network` now holds an instance of the `ContainedEntityNetwork`,
which we can also pass to the `GmshWriter`:

```cpp
GmshWriter writer(network);
writer.write("network", // filename of the .geo files (will add extension .geo automatically)
             0.1,       // element size to be used on the quadrilaterals
             0.2);      // element size to be used in the domain away from the quadrilaterals
```

As you can see, one can specify different mesh sizes to be used on the fracture
entities and in the rest of the domain.

[0]: https://git.iws.uni-stuttgart.de/DennisGlaeser/frackit/tree/master/appl/example1
[1]: https://git.iws.uni-stuttgart.de/DennisGlaeser/frackit/tree/master/geometry/cylinder.hh
[2]: https://git.iws.uni-stuttgart.de/DennisGlaeser/frackit/tree/master/appl/example3
