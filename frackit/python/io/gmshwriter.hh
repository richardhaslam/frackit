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
#ifndef FRACKIT_PYTHON_IO_GMSH_WRITER_HH
#define FRACKIT_PYTHON_IO_GMSH_WRITER_HH

#include <frackit/io/gmshwriter.hh>
#include <frackit/entitynetwork/entitynetwork.hh>
#include <frackit/entitynetwork/containedentitynetwork.hh>

namespace Frackit::Python {

namespace py = pybind11;

template<class ctype>
void registerGmshWriter(py::module& module)
{
    using EntityNetwork = Frackit::EntityNetwork;
    using ContainedEntityNetwork = Frackit::ContainedEntityNetwork;

    py::class_<GmshWriter> cls(module, "GmshWriter");
    cls.def(py::init<const EntityNetwork&>());
    cls.def(py::init<const ContainedEntityNetwork&>());

    using namespace py::literals;
    cls.def("write",
            py::overload_cast<const std::string&, ctype, ctype>(&GmshWriter::template write<ctype>),
            "fileName"_a, "meshSizeAtEntities"_a, "meshSizeAtBoundary"_a,
            "write the entity network to a .geo file");
    cls.def("write",
            py::overload_cast<const std::string&, ctype>(&GmshWriter::template write<ctype>),
            "fileName"_a, "meshSizeAtEntities"_a,
            "write the entity network to a .geo file");
}

} // end namespace Frackit::Python

#endif
