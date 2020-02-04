---
title: 'Frackit: a framework for stochastic fracture network generation and analysis'
tags:
  - Fractures
  - Fracture network
  - Porous media
authors:
  - name: Dennis Gläser
    orcid: 0000-0001-9646-881X
    affiliation: 1
  - name: Bernd Flemisch
    affiliation: 1
  - name: Holger Class
    affiliation: 1
  - name: Rainer Helmig
    affiliation: 1
affiliations:
 - name: Department of Hydromechanics and Modelling of Hydrosystems, University of Stuttgart, Pfaffenwaldring 61, 70569 Stuttgart, Germany
   index: 1
date: 30 January 2020
bibliography: paper.bib
---


# Summary

The numerical simulation of flow and transport phenomena in fractured porous media
is an active field of research, given the importance of fractures in many
geotechnical engineering applications, as for example groundwater management
[@qian2014], enhanced oil recovery techniques [@torabi2012],
geothermal energy [@mcfarland1976; @shaik2011] or unconventional
natural gas production [@sovacool2014]. A number of mathematical models and numerical
schemes, aiming at an accurate description of flow through fractured rock, have been
presented recently
(see e.g. @ahmed2015; @ahmed2017; @brenner2018; @koppel2019; @schadle2019; @nordbotten2019). Many of these describe the fractures as lower-dimensional
geometries, that is, as curves or planes embedded in two- or three-dimensional
space, respectively. On those, integrated balance equations are solved together
with transmission conditions describing the interaction with the surrounding medium.
Moreover, it is often required that the computational meshes used for the different
domains are conforming in the sense that the faces of the discretization used for
the bulk medium coincide with the discretization of the fractures (see image below).

![Exemplary grids used with numerical schemes that require conformity of the discretizations.](doc/img/examplegrids.png)

Information on the in-situ locations of fractures is typically sparse and
difficult to determine. In response to this, a common approach is to study
the hydraulic properties of rock in function of the fracture network topology by means
of numerical simulations performed on stochastically generated fracture networks.
Such investigations have been presented, among others, in
@Kazumasa2003; @Assteerawatt2008; @lee2015fracture; @zhang2015finite; @lee2019stochastic.
An open-source Matlab code for the stochastic generation and analysis fracture
networks in two- and three-dimensional space has been presented in @alghalandis2017.
However, the code is limited to linear (polygonal) fracture geometries, embedded
in hexahedral domains.

``Frackit`` is a C++-framework for the stochastic generation of fracture networks
composed of polygonal and/or elliptical geometries, embedded in arbitrary domain
shapes. It makes extensive use of the open-source Computed-Aided-Design (CAD)
library [OpenCascade][4] ([opencascade.com][5]), which offers great flexibility with respect to the
geometries that can be used. Moreover, a large number of standard CAD file formats
is supported for input/output of geometrical shapes. This allows users of ``Frackit``
to read in externally generated domain shapes
(for instance, from measurements and/or created using CAD software),
and to generate fracture networks within these domains. Output routines to standard
file formats enable users to then construct computational meshes of the generated
geometries using a variety of tools. In particular, ``Frackit`` offers
output routines to the (.geo) file format used by [Gmsh][2] [@gmsh2009],
which is an open-source mesh generator that is widely used in academic
research (see e.g. @keilegavlen2017; @berge2019).

The geometric data produced by ``Frackit`` contains the complete fragmentation
of all geometric entities involved, i.e. the intersection geometries between
all entities are computed. Thus, this information can be directly used in the
context of discrete fracture-matrix (dfm) simulations in a conforming way as
described above. For instance, the open-source simulator [DuMuX][3] [@Dumux; @koch2019dumux31]
contains a module for conforming dfm simulations of single- and multi-phase
flow through fractured porous media, which has been used in several works
[@glaeser2017, @glaeser2019; @andrianov2019].
It  supports the [Gmsh][2] file format (.msh), and thus, ``Frackit`` can be used in
a fully open-source toolchain with [Gmsh][2] and [DuMuX][3] to generate random
fracture networks, construct computational meshes, and perform analyses on them
by means of numerical simulations.

The design of ``Frackit`` is such that there is no predefined program flow, but
instead, users should implement their own applications using the provided classes
and functions, which allows for full customization of each step of the network
generation. Besides this, in the case of available measurement data, one could
skip the network generation process and use ``Frackit`` to compute the fragmentation
of the measured data and to generate CAD files for subsequent meshing.

# Concept

The functionality provided in ``Frackit`` follows from a conceptual division of
the network generation into three basic steps:

* Random generation of raw fracture entities based on statistical parameters
* Evaluation of geometric constraints for a new entity candidate against
  previously generated entities
* Fragmentation of the generated raw entities and the embedding domain

The last of these steps and the motivation for it has been discussed
above. In the following, we want to discuss the other two steps in more detail.

## Random generation of raw fracture entities

In the network generation procedure, a domain is populated with fracture entities
that are generated following user-defined statistical properties regarding their
size, orientation and spatial distribution. In ``Frackit``, this process is termed
_geometry sampling_ and is realized in the code in _sampler_ classes. In the current
implementation, there are two such sampler classes available, which sample quadrilaterals
and elliptical disks in three-dimensional space. A sampler class of ``Frackit``
receives an instance of a `PointSampler` implementation and a number of probability
distributions that define the size and orientation of the raw entities.
`PointSampler` classes are used to sample the spatial distribution of the geometries
inside a domain geometry. For example, a point sampler that samples points uniformly
within the unit cube (defined in the variable `domain`) could be constructed like this:

```cpp
// the type used for coordinates values
using ctype = double;

// define axis-aligned box in which to sample the centers points
using Domain = Frackit::Box<ctype>;
Domain domain(0.0, 0.0, 0.0,  // xmin, ymin, zmin
              1.0, 1.0, 1.0); // xmax, ymax, zmax

// let us uniformly sample points within this box
const auto pointSampler = Frackit::makeUniformPointSampler(domain);
```

The convenience function `makeUniformPointSampler()` can be used for uniform
sampling over the provided domain geometry. For nun-uniform samplers, one can write

```cpp
const auto pointSampler = Frackit::makePointSampler<Traits>(domain);
```

where in the `Traits` class users define the type of distribution to be used
for each coordinate direction.
Inside a geometry sampler class, a geometry is created by sampling a
point from the point sampler, and then constructing a geometry around this
point using the provided distributions for its size and orientation. For example,
the `QuadrilateralSampler` class expects distributions for the strike angle,
dip angle, edge length and a threshold value for the minimum allowed edge length.
The following piece of code shows how an instance of the `QuadrilateralSampler`
class, using uniform distributions for all parameters regarding orientation and
size, can be created (we reuse the `pointSampler` variable defined in the previous
code snippet):

```cpp
// let us use uniform distributions for the quadrilateral parameters
using Distro = std::normal_distribution<ctype>;

// Distributions for strike & dip angle & edge length
Distro strikeAngleDistro(toRadians(45.0), // mean value
                         toRadians(5.0)); // standard deviation
Distro dipAngleDistro(toRadians(45.0), // mean value
                     toRadians(5.0));  // standard deviation
Distro edgeLengthDistro(0.5,  // mean value
                        0.1); // standard deviation

// instance of the quadrilateral sampler class
using QuadSampler = Frackit::QuadrilateralSampler</*spaceDimension*/3>;
QuadSampler quadSampler(pointSampler,
                        strikeAngleDistro,
                        dipAngleDistro,
                        edgeLengthDistro,
                        0.05); // threshold for minimum edge length
```

As for point samplers, one can use non-uniform distributions by implementing
a `Traits` class which is then passed to the `QuadrilateralSampler` as template
argument. The definitions of the strike and dip angles as used within the
`QuadrilateralSampler` class are illustrated in the figure below. Consider a
quadrilateral whose center is the origin and which lies in the plane defined by
the two basis vectors $\mathbf{b}_1$ and $\mathbf{b}_2$. The latter lies in the
$x$-$y$-plane and the strike angle is the angle between the $y$-axis and $\mathbf{b}_2$.
The dip angle describes the angle between $\mathbf{b}_1$ and the $x$-$y$-plane.

![Illustration of the strike and dip angles involved in the random generation of quadrilaterals. The grey plane with the structured mesh illustrates the $x$-$y$-plane.](doc/img/quadsampler.png)

In the code, random generation of geometries from sampler classes occurs by using
the `()` operator. For example, from the `quadSampler` variable defined in the
previous code snippet, we obtain a random quadrilateral by writing:

```cpp
// generate random quadrilateral
const auto quad = quadSampler();
```

## Evaluation of geometric constraints

While the domain is populated with the raw fracture entities, users have the
possibility to enforce geometric constraints between different entities in order
to enforce topological characteristics as e.g. fracture spacing. Besides this,
constraints can be used to avoid very small length scales that could cause problems
during mesh generation or could lead to ill-shaped elements. In the code, constraints
can be defined and evaluated using the `EntityNetworkConstraints` class. These have
to be fulfilled by a new fracture entity candidate against previously accepted
entities. If any of the defined constraints is violated, the candidate may be
rejected and a new one is sampled. The current implementation of the
`EntityNetworkConstraints` class allows users
to define a minimum distance between two entities that do not intersect. If two
entities intersect, one can choose to enforce a minimum length of the intersection
curve, a minimum intersection angle and a minimum distance between the intersection
curve and the boundaries of the intersecting entities. An illustration of this
is shown in the figure below.

![Overview of the geometric constraints that can be defined between entities.](doc/img/constraints.png)

The following code snippet illustrates how to set up an instance of the
`EntityNetworkConstraints` class:

```cpp
// Instantiate constraints class. This leaves all constraints deactivated.
Frackit::EntityNetworkConstraints<ctype> constraints;

// Set values to activate constraints
constraints.setMinDistance(0.1);               // in meter
constraints.setMinIntersectingAngle(M_PI/4.0); // in radians
constraints.setMinIntersectionMagnitude(0.05); // in meter
constraints.setMinIntersectionDistance(0.05);  // in meter

// define tolerance value to be used for intersections
constraints.setIntersectionEpsilon(1e-6);
```

When using the default constructor of `EntityNetworkConstraints`, all constraints
are inactive, and when defining values for the different constraint types, these
get activated internally. Moreover, one can define the tolerance value that should
be used in the intersection algorithms between entities. If no tolerance is set,
a default tolerance is computed based on the size of the entities for which the
intersection is to be determined. For two quadrilaterals `quad1` and `quad2`, one
can then evaluate the defined constraints by writing:

```cpp
bool fulfilled = constraints.evaluate(quad1, quad2);
```

The function `evaluate()` returns true if all constraints are fulfilled. One can
also check the fulfillment of the constraints of a new candidate against an entire
set of entities. Let `quad` be a new candidate for a quadrilateral, and
`quadSet` be a vector of quadrilaterals (`std::vector< Quadrilateral<ctype> >`),
then one can write

```cpp
bool fulfilled = constraints.evaluate(quadSet, quad);
```

to evaluate the constraints between `quad` and all entities stored in `quadSet`.


# Example application

In this exemplary application we want to briefly outline what the workflow using
`Frackit` together with [Gmsh][2] and [DuMuX][3] could look like. The images are
taken from the `Frackit` documentation ([git.iws.uni-stuttgart.de/tools/frackit][0])
and the configurations of the geometry samplers are, apart from small modifications,
very similar to the ones used in [example 3][1] provided in the `Frackit` repository.
For further details on how to set up such configurations we refer to the source code
and the documentation of that example in the repository.

Let us consider a domain consisting of three solid layers, of which we want to
generate a fracture network only in the center volume. With the following piece
of code we read in the domain geometry from the provided file, extract the three
volumes of it and select the middle one as the one in which we want to place the
fracture network.

```cpp
/////////////////////////////////////////////////////
// 1. Read in the domain geometry from .brep file. //
//    The file name is defined in CMakeLists.txt   //
/////////////////////////////////////////////////////
const auto domainShape = Frackit::OCCUtilities::readShape(BREPFILE);

// obtain the three solids contained in the file
const auto solids = Frackit::OCCUtilities::getSolids(domainShape);

// The sub-domain we want to create a network in is the center one.
const auto& networkDomain = solids[1];

// get the bounding box of the domain
const auto bBox = Frackit::OCCUtilities::getBoundingBox(networkDomain);
```

The last command constructs the bounding box of the center volume of our domain,
which we can then use to instantiate point sampler classes that define
the spatial distribution of the fracture entities. With these, we can construct
geometry samplers as outlined above. In this example, we define three
geometry sampler instances to sample from three different orientations of fractures,
and we use quadrilaterals for two of the orientations and elliptical disks for
the third orientation. Moreover, we define different constraints that should be
fulfilled between the entities of different orientations. As mentioned above, details
on how to implement such settings can be found in [example 3][1] in the `Frackit` repository.

A number of fractures is then generated for each orientation. Subsequently, the
raw entities and the three volumes of the domain are cast into an instance
of the `ContainedEntityNetwork` class. This can be used to define arbitrarily many
(sub-)domains, and to insert entities to be embedded in a specific sub-domain.
The `ContainedEntityNetwork` computes and stores the fragments of all entities
and sub-domains resulting from mutual intersection. Output routines for instances
of this class are implemented, which generate geometry files that are ready to be
meshed using designated tools, for example, [Gmsh][2].

The image below illustrates the workflow chosen in this example, using ``Frackit``
to generate a random fracture network, [Gmsh][2] to mesh the resulting geometry,
and [DuMuX][3] to perform a single-phase flow simulation on the resulting mesh.
The bottom picture shows the pressure distribution on the fractures and the
velocities in the domain as computed with [DuMuX][3], using the illustrated
boundary conditions.

![Illustration of the workflow using Frackit, Gmsh and DuMuX in the exemplary application.](doc/img/network_and_solution.png)

The source code of this example, including installation instructions, can be found
at https://git.iws.uni-stuttgart.de/dumux-pub/glaeser2020a.

# Future developments

We are planning to add fracture network characterization capabilities, such as
the detection of isolated clusters of fractures or the determination of connectivity
measures. In order to do this efficiently, we want to integrate data
structures and algorithms for graphs, together with functionalities to translate
the generated fracture networks into graph representations. Besides this, we want
to develop python bindings for ``Frackit`` to provide an easy-to-use and high-level
interface.


# Acknowledgements

We thank the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation)
for supporting this work by funding SFB 1313, Project Number 327154368.

# References


[0]: https://git.iws.uni-stuttgart.de/tools/frackit
[1]: https://git.iws.uni-stuttgart.de/tools/frackit/tree/master/appl/example3
[2]: http://www.gmsh.info/
[3]: https://dumux.org/
[4]: https://www.opencascade.com/content/download-center
[5]: https://www.opencascade.com
