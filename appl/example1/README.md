<!--- Example picture --->
<p align="center">
    <img src="../../doc/img/example1_network.png" alt="frackit example 1" width="800"/>
</p>

Example 1
=========

In this exemplary application, a network of quadrilateral fractures with two main
orientations is generated within the unit cube (see image above). The main file
containing the source code to this example is the file `example1.cc` which is
located in this folder. Note that this description focuses on the C++ implementation,
but in `example1.py` you can find how to realize this example using the Frackit
python bindings.

<b> _In this example, you will learn how to:_ </b>

* use the sampler class for quadrilaterals
* define and enforce geometric constraints between entities
* construct an __EntityNetwork__ out of the raw entities to gather connectivity information


### Quadrilateral samplers

Two main orientations are considered for the quadrilaterals, and for both
a corresponding instance of the `QuadrilateralSampler` class is created.
For example, we instantiate an instance of this class by writing:

```cpp
static constexpr int worldDimension = 3;

using ctype = double;
using NormalDistro = std::normal_distribution<ctype>;
using UniformDistro = std::uniform_real_distribution<ctype>;
using QuadSampler = QuadrilateralSampler<worldDimension>;

Box<ctype> domain(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
QuadSampler quadSampler(makeUniformPointSampler(domain),               // point sampler that samples the center points of the quadrilaterals
                        NormalDistro(toRadians(0.0), toRadians(5.0)),  // strike angle: mean value & standard deviation
                        NormalDistro(toRadians(0.0), toRadians(5.0)),  // dip angle: mean value & standard deviation
                        UniformDistro(0.4, 0.8),                       // strike length
                        UniformDistro(0.4, 0.8));                      // dip length
```

The first constructor argument is a point sampler with which the center points of
the quadrilaterals are sampled. Here we use uniformly sampled points in the unit
cube, which is represented by an instance of the `Box` class,  stored in the
variable `domain`. The second and third arguments define the distributions for
the strike and dip angles (for details see the [class documentation][2]), where in
this case, we use uniform distributions with
a mean value of 0° and a standard deviation of 5°. The last two arguments
are distributions for the sizes of the quadrilaterals in strike and dip direction.

In the example, the quadrilaterals are sampled from the two samplers `quadSampler1` and
`quadSampler2`, using the `()` operator:

```cpp
auto quad = sampleIntoSet1 ? quadSampler1() : quadSampler2();
```

Here, we use the boolean variable `sampleIntoSet1` to determine from
which sampler we should sample the next quadrilateral (more details follow below).
The variable `quad` holds a new candidate for an entity of the network, however,
we want to enforce certain constraints such as a minimum distance between
entities.

### Constraints definitions

In order for certain constraints to be fulfilled among the entities, we use instances
of the `EntityNetworkConstraints` class and
configure it as desired. For example, the constraints on entities of the same
orientation are defined in this example as follows:

```cpp
// We want to enforce some constraints on the set of quadrilaterals.
// In particular, for entities of the same set we want a minimum spacing
// distance of 5cm, and the quadrilaterals must not intersect in angles
// less than 30°. Moreover, if they intersect, we don't want intersection
// edges whose length is smaller than 5cm, and, the intersection should not
// be too close to the boundary of one of two intersecting quadrilaterals. Here: 5cm.
EntityNetworkConstraints<ctype> constraintsOnSelf;
constraintsOnSelf.setMinDistance(0.05);
constraintsOnSelf.setMinIntersectingAngle(toRadians(30.0));
constraintsOnSelf.setMinIntersectionMagnitude(0.05);
constraintsOnSelf.setMinIntersectionDistance(0.05);
```

Among entities of different orientation, we want to ensure larger intersection angles.
To this end, we simply create a copy of the constraints object, but define a larger
minimum intersection angle:

```cpp
// with respect to entities of the other set, we want to have larger intersection angles
 auto constraintsOnOther = constraintsOnSelf;
constraintsOnOther.setMinIntersectingAngle(toRadians(40.0));
```

### Entity sampling

In the main loop of quadrilateral generation, we sample candidates for new
entities as described above. Then, the fulfilment of the above-mentioned
constraints is evaluated against the other quadrilaterals with:

```cpp
auto& entitySet = sampleIntoSet1 ? entitySet1 : entitySet2;
if (!constraintsOnSelf.evaluate(entitySet, quad))
{ status.increaseRejectedCounter(); continue; }
```

where `entityset1` and `entitySet2` are of type `std::vector<Quadrilateral>` and
store all quadrilaterals that are accepted. The function `evaluate` of the
`EntityNetworkConstraints` class evaluates the constraints for `quad` against all
entities contained in `entitySet` and returns `true` only if no violation of
any of the defined constraints has been detected. After an admissible quadrilateral
has been generated, we add it to the current set by

```cpp
// the quadrilateral is admissible
entitySet.push_back(quad);
```

and make sure that in the next iteration we sample a candidate of the other orientation:

```cpp
// sample into the other set the next time
sampleIntoSet1 = !sampleIntoSet1;
```

In [Example 3][0] we will get to know how to use helper classes
that store different entity sets and automatically sample from various sampler
classes such that this can be written more easily.

### Network construction

After the desired number of entities has been generated, the entities are cast
into an entity network using the builder class:

```cpp
EntityNetworkBuilder<ctype> builder;
builder.addEntities(entitySet1);
builder.addEntities(entitySet2);

const auto network = builder.build();
```

This network can then be written to disk, for example in [Gmsh][1] (.geo) file format:

```cpp
GmshWriter writer(network);
writer.write("network", // filename of the .geo files (will add extension .geo automatically)
             0.1);      // element size to be used on the network entities
```

Note that with the `EntityNetworkBuilder` class we have created a network that
solely carries information about the fracture entities. We have not defined any
domain in this example, the unit cube in the variable `domain` was only used to
sample the center points of the quadrilaterals. Thus, the geometry files written
by the `GmshWriter` also only contain data on the fracture entities. This can be
used in contexts where one is only interested in the fractures. In the following
examples we will see how to construct fracture networks embedded in one or more
(sub-)domains.

| [:arrow_right: Go to example 2](https://git.iws.uni-stuttgart.de/tools/frackit/tree/master/appl/example2) |
|---:|

[0]: https://git.iws.uni-stuttgart.de/tools/frackit/tree/master/appl/example3/README.md
[1]: http://gmsh.info/
[2]: https://git.iws.uni-stuttgart.de/tools/frackit/blob/master/frackit/sampling/quadrilateralsampler.hh
