// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef FRACKIT_PYTHON_ENTITYNETWORK_BUILDER_HH
#define FRACKIT_PYTHON_ENTITYNETWORK_BUILDER_HH

#include <pybind11/pybind11.h>

// (currently) supported domain geometries
#include <frackit/geometry/box.hh>
#include <frackit/geometry/cylinder.hh>

// (currently) supported entity geometries
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/polygon.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/python/geometry/brepwrapper.hh>

#include <frackit/common/id.hh>
#include <frackit/entitynetwork/networkbuilder.hh>
#include "networkbuilderwrappers.hh"

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    template<class Domain, class Builder, class... Bases>
    void registerSubDomainAdder(py::class_<Builder, Bases...>& cls)
    {
        cls.def("addSubDomain", &Builder::template addSubDomain<Domain>, "add subdomain with the given id");
        cls.def("addConfiningSubDomain", &Builder::template addConfiningSubDomain<Domain>, "add confining subdomain with the given id");
    }

    template<class Entity, class Builder, class... Bases>
    void registerEntityAdder(py::class_<Builder, Bases...>& cls)
    {
        cls.def("addEntity",
                &Builder::template addEntity<Entity>,
                "add entity to the network");
    }

    template<class Entity, class Builder, class... Bases>
    void registerSubDomainEntityAdder(py::class_<Builder, Bases...>& cls)
    {
        cls.def("addSubDomainEntity",
                &Builder::template addSubDomainEntity<Entity>,
                "add entity to be embedded in the subdomain with the given id");
    }

    template<class ctype, class Builder, class... Bases>
    void registerSubDomainAdders(py::class_<Builder, Bases...>& cls)
    {
        registerSubDomainAdder< Box<ctype> >(cls);
        registerSubDomainAdder< Cylinder<ctype> >(cls);
        registerSubDomainAdder< SolidWrapper >(cls);
    }

    template<class ctype, class Wrapper, class... Bases>
    void registerSubDomainEntityAdders(py::class_<Wrapper, Bases...>& cls)
    {
        registerSubDomainEntityAdder< Disk<ctype> >(cls);
        registerSubDomainEntityAdder< Quadrilateral<ctype, 3> >(cls);
        registerSubDomainEntityAdder< Polygon<ctype, 3> >(cls);
        registerSubDomainEntityAdder< FaceWrapper >(cls);
    }

    template<class ctype, class Wrapper, class... Bases>
    void registerEntityAdders(py::class_<Wrapper, Bases...>& cls)
    {
        registerEntityAdder< Disk<ctype> >(cls);
        registerEntityAdder< Quadrilateral<ctype, 3> >(cls);
        registerEntityAdder< Polygon<ctype, 3> >(cls);
        registerEntityAdder< FaceWrapper >(cls);
    }

} // end namespace detail

template<class ctype>
void registerEntityNetworkBuilders(py::module& module)
{
    // register common base class
    using Base = EntityNetworkBuilderBase<ctype>;
    py::class_<Base> base(module, "_EntityNetworkBuilderBase");
    base.def("setConfineToSubDomainUnion",
             &Base::setConfineToSubDomainUnion,
             "set if entities are to be confined by the union of all subdomains (for non-confining subdomains)");

    base.def("setEpsilon", &Base::setEpsilon, "define tolerance value for intersections");
    base.def("setDefaultEpsilon", &Base::setDefaultEpsilon, "restore default tolerance value for intersections");
    base.def("clear", &Base::clear, "clears all inserted data");

    // register base classes for builders
    using Builder = EntityNetworkBuilder<ctype>;
    py::class_<Builder, Base> builder(module, "_EntityNetworkBuilder");
    builder.def("build", &Builder::build, "build the network from the inserted entities");

    using ContainedBuilder = ContainedEntityNetworkBuilder<ctype>;
    py::class_<ContainedBuilder, Base> containedBuilder(module, "_ContainedEntityNetworkBuilder");
    containedBuilder.def("build", &ContainedBuilder::build, "build the contained network from the inserted entities");

    // register wrapper classes
    using BuilderWrapper = Detail::EntityNetworkBuilderWrapper<ctype>;
    py::class_<BuilderWrapper, Builder> builderWrapper(module, "_EntityNetworkBuilderWrapper");
    builderWrapper.def(py::init());
    Detail::registerSubDomainEntityAdders<ctype>(builderWrapper);
    Detail::registerEntityAdders<ctype>(builderWrapper);

    // register adder functions for (sub-)domains
    Detail::registerSubDomainAdders<ctype>(builderWrapper);

    using ContainedBuilderWrapper = Detail::ContainedEntityNetworkBuilderWrapper<ctype>;
    py::class_<ContainedBuilderWrapper, ContainedBuilder> containedBuilderWrapper(module, "_ContainedEntityNetworkBuilderWrapper");
    containedBuilderWrapper.def(py::init());
    Detail::registerSubDomainEntityAdders<ctype>(containedBuilderWrapper);

    // register adder functions for (sub-)domains
    Detail::registerSubDomainAdders<ctype>(containedBuilderWrapper);
}

} // end namespace Frackit::Python

#endif
