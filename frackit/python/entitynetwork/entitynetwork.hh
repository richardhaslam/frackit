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
#ifndef FRACKIT_PYTHON_ENTITY_NETWORK_HH
#define FRACKIT_PYTHON_ENTITY_NETWORK_HH

#include <pybind11/pybind11.h>

#include <frackit/entitynetwork/entitynetwork.hh>
#include <frackit/entitynetwork/containedentitynetwork.hh>

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    void registerEntityNetwork(py::module& module)
    {
        py::class_<EntityNetwork>(module, "EntityNetwork");
    }

    void registerContainedEntityNetwork(py::module& module)
    {
        py::class_<ContainedEntityNetwork>(module, "ContainedEntityNetwork");
    }

} // end namespace Detail

void registerEntityNetworks(py::module& module)
{
    Detail::registerEntityNetwork(module);
    Detail::registerContainedEntityNetwork(module);
}

} // end namespace Frackit::Python

#endif
