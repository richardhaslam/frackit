<!--- Example picture --->
<p align="center">
    <img src="doc/img/titlepicture.png" alt="frackit example" width="800"/>
</p>

What is Frackit?
================

Frackit is an open-source framework for the stochastic generation of fracture networks.
The implementation is written in C++ language and extensively uses
[OpenCascade][2], an open-source Computer Aided Design (CAD) library.

The geometries generated by Frackit can be exported into a number of standard CAD
file formats supported by [OpenCascade][2]. This allows for visualization with a
variety of tools, as for example the [CAD Assistant][3] developed by [OpenCascade][2],
or [Gmsh][1], a three-dimensional finite element mesh generator.

Coupling to [DuMuX][0]
======================

Frackit also features output routines to [Gmsh][1] file format (.geo), where mesh
size specifications can be defined upon fracture network generation. This allows
for the generation of computational meshes using [Gmsh][1], which are supported
by the open-source simulator [DuMuX][0] for flow and transport in porous media.
The .geo files produced by Frackit lead to three-dimensional meshes that are aligned
with the fracture geometries such that the element faces coincide with the fractures,
which can be directly plugged into the [DuMuX][0] module for discrete fracture-matrix
simulations of single- and two-phase flow (see e.g. https://arxiv.org/pdf/1909.05052.pdf).

For example, the pressure distribution on the fractures and the fluid velocities
in the domain, for the fracture network shown above, could look like this:

<p align="center">
    <img src="doc/img/pressuresolution.png" alt="dumux frackit example" width="800"/>
</p>

The source code to this example, including installation instructions, can be found
[here][13].

General Concept
===============

### Geometry Sampling
The generation of the fracture networks occurs by randomly sampling instances of the
desired fracture geometry on the basis of probability distribution functions, which
can be defined by the user. The implementation allows for both selecting the type of
distribution (uniform, exponential, etc.) as well as the distribution parameters.

### Constraints Evaluation
After the generation of a new candidate for a fracture entity, a number of constraints can
be evaluated for it. These can be used to enforce topological characteristics of the
fracture network, e.g. fracture spacing, by defining a minimum distance between entities.
Other constraints are targeted mainly at guaranteeing certain mesh properties by avoiding
very small length scales, and thus, small elements in the computational mesh. As constraints
one can define a minimum length scale of the intersections between fracture entities, a
minimum intersection angle, and a minimum distance of the intersection geometry to the boundary
of the intersecting entities. If the user-defined constraints are not fulfilled, the candidate
is rejected.

### Fragmentation of the network
After the desired number of fracture entities have been generated, an __EntityNetwork__
can be constructed from the raw entities. This intersects and fragments all entities, and
if desired, one can confine the network to a domain of choice.

Note that none of these steps are not mandatory, nor do they occur in a fixed order.
Instead, users are motivated to implement their own applications using the provided
functions, classes and concepts. The modular design of the above-mentioned building
blocks enables users to fully customize each step to their needs. Three exemplary
applications are contained in this repository:

* [Example 1][5] Generation of a simple network consisting of quadrilaterals with two main orientations.
* [Example 2][6] Generation of a network of quadrilaterals embedded in a cylindrical domain (mimicking a core sample).
* [Example 3][4] Generation of a network consisting of both disks and quadrilaterals, confined and contained in one layer of a domain that is composed of three layers.


Documentation
=============

A class documentation can be generated from the source code using
[Doxygen][8] (see Installation notes).
Moreover, the [Examples][7] contained in this repository provide a good overview over
the capabilities of Frackit and can serve as a starting point to develop your own
application.


Installation
============

Please note that the following requirements need to be installed:

* OpenCascade (>= 7.3.0)
* CMake (>3.0)
* C, C++ compiler (C++17 required)
* Optional: Doxygen (>= 1.8)

### Installation of OpenCascade
Frackit requires the [OpenCascade][2] library to be installed on your system.
You can download the source code [HERE][2],
and details on the installation can be found [HERE][10].
Please note that [OpenCascade][2] requires further 3rd party libraries, of which
the mandatory ones are Tcl/Tk and FreeType (see this [link][11]).
On Ubuntu, both of these can be installed from the command line.

### Building Frackit under Linux
After [OpenCascade][2] and the other requirements listed above have been installed,
clone this repository within your folder of choice by typing:

```sh
git clone https://git.iws.uni-stuttgart.de/DennisGlaeser/frackit.git
```

Then, create the build directory in which you want to compile the applications
and run cmake from it. For instance:

```sh
mkdir build
cd build
cmake ../
```

If cmake cannot find your installation of [OpenCascade][2], you probably installed it
into a non-standard location. In this case, you can define __HINTS__ for cmake to search
for it. In particular, you would have to change the line

```cmake
find_path(OCC_INC "Standard_Version.hxx" PATH_SUFFIXES opencascade include/opencascade)
```

to

```cmake
find_path(OCC_INC "Standard_Version.hxx" HINTS MY_OCC_INCLUDE_FOLDER)
```

in the _CMakeLists.txt_ file of the top folder of Frackit, substituting
MY_OCC_INCLUDE_FOLDER with the path to the source files of [OpenCascade][2]
on your system. The same has to be done for the required packages of
[OpenCascade][2], i.e. in the line

```cmake
find_library(OCC_LIB ${OCC}
             PATH_SUFFIXES lib ${OCC_SYS_NAME}/lib ${OCC_SYS_NAME}/vc8/lib)
```

you can define HINTS for cmake to find your installation folder of [OpenCascade][2].
Once cmake finished successfully, you could now compile the class documentation:

```sh
make doc_doxygen
```

and open it using a web browser, for example _chrome_:

```sh
google-chrome doc/doxygen/html/index.html
```

To build and run all unit tests, type

```sh
make build_tests
ctest
```

Moreover, you can build all example applications by typing

```sh
make build_example_applications
```

If this compiled successfully, you can go to the folder of the first example and
run it by typing:

```sh
cd appl/example1
./example1
```


License
=======

Frackit is licensed under the terms and conditions of the GNU General
Public License (GPL) version 3 or - at your option - any later
version. The GPL can be [read online][9] or in the [LICENSE.md](LICENSE.md) file
provided in the top folder of Frackit. See the file [LICENSE.md](LICENSE.md) for
full copying permissions.


Contributing
=============

Contributions are highly appreciated.
For bug reports, please file an [issue](https://git.iws.uni-stuttgart.de/DennisGlaeser/frackit/issues).
If you want to contribute with new features of improvements of existing code, please

* Get an account for our GitLab instance at https://git.iws.uni-stuttgart.de/
* Fork this project
* Push your changes to your fork on some branch.
* Open a merge request using the branch you pushed your changes to as the source branch and the master
of the Frackit repository as the target branch.
* Follow the discussion on the merge request to see what improvements should be done to the branch
before merging.
* If you have developer status you don't need to do a fork and you can create branches directly.

In this project, we follow the [styleguide][2] of the [DuMuX][0] project.
Please have a look at these before you start coding your contributions.


Ongoing Developments
====================

* Kernel for fracture network characterization (e.g. connectivity)
* Capabilities for translation of fracture networks into graphs, allowing for efficient characterization of the network based on its graph representation.
* Python bindings

[0]: https://dumux.org
[1]: http://gmsh.info/
[2]: https://www.opencascade.com/content/download-center
[3]: https://www.opencascade.com/content/cad-assistant
[4]: https://git.iws.uni-stuttgart.de/DennisGlaeser/frackit/tree/master/appl/example3
[5]: https://git.iws.uni-stuttgart.de/DennisGlaeser/frackit/tree/master/appl/example1
[6]: https://git.iws.uni-stuttgart.de/DennisGlaeser/frackit/tree/master/appl/example2
[7]: https://git.iws.uni-stuttgart.de/DennisGlaeser/frackit/tree/master/appl/
[8]: http://www.doxygen.org/index.html
[9]: https://www.gnu.org/licenses/gpl-3.0.en.html

[10]: https://www.opencascade.com/doc/occt-6.9.1/overview/html/occt_dev_guides__building_cmake.html
[11]: https://www.opencascade.com/doc/occt-6.9.1/overview/html/occt_dev_guides__building_3rdparty_linux.html
[12]: https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/doc/styleguide.md
[12]: TODO:DUMUXPUB
