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
#ifndef FRACKIT_PYTHON_ENTITYNETWORK_CONSTRAINTS_HH
#define FRACKIT_PYTHON_ENTITYNETWORK_CONSTRAINTS_HH

#include <algorithm>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <frackit/geometry/disk.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/entitynetwork/constraints.hh>

#include <frackit/python/occutilities/brepwrapper.hh>

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    // Wrapper class to resolve overload issues for evaluate()
    template<class ST>
    class EntityNetworkConstraintsWrapper
    : public EntityNetworkConstraints<ST>
    {
        using ParentType = EntityNetworkConstraints<ST>;

    public:
        // evaluate the constraint between two geometries
        template<class Geo1, class Geo2>
        bool evaluateBinary(const Geo1& geo1, const Geo2& geo2) const
        {
            static constexpr auto isWrapper1 = OCCUtilities::IsBRepWrapper<Geo1>::value;
            static constexpr auto isWrapper2 = OCCUtilities::IsBRepWrapper<Geo2>::value;

            if constexpr (!isWrapper1 && !isWrapper2)
                return ParentType::evaluate(geo1, geo2);
            else if constexpr (!isWrapper1 && isWrapper2)
                return ParentType::evaluate(geo1, geo2.get());
            else if constexpr (isWrapper1 && !isWrapper2)
                return ParentType::evaluate(geo1.get(), geo2);
            else
                return ParentType::evaluate(geo1.get(), geo2.get());
        }
    };

    template<class Constraints>
    void registerSetters(py::class_<Constraints>& cls)
    {
        cls.def("setMinDistance", &Constraints::setMinDistance,
                "define the minimum distance between entities");
        cls.def("setMinIntersectingAngle", &Constraints::setMinIntersectingAngle,
                "define the minimum angle between intersecting entities");
        cls.def("setMinIntersectionMagnitude", &Constraints::setMinIntersectionMagnitude,
                "define the minimum magnitude of intersections between entities");
        cls.def("setMinIntersectionDistance", &Constraints::setMinIntersectionDistance,
                "define the minimum distance between an intersection and the boundaries of the intersecting entities");
        cls.def("neglectMinDistance", &Constraints::neglectMinDistance,
                "deactivate minimum distance constraint");
        cls.def("neglectMinIntersectionAngle", &Constraints::neglectMinIntersectionAngle,
                "deactivate minimum intersection angle constraint");
        cls.def("neglectMinIntersectionMagnitude", &Constraints::neglectMinIntersectionMagnitude,
                "deactivate minimum intersection magnitude constraint");
        cls.def("neglectMinIntersectionDistance", &Constraints::neglectMinIntersectionDistance,
                "deactivate constraint on minimum distance between intersections and entity boundaries");
        cls.def("setIntersectionEpsilon", &Constraints::setIntersectionEpsilon,
                "define a tolerance value to be used in intersection computations");
        cls.def("setDefaultIntersectionEpsilon", &Constraints::setDefaultIntersectionEpsilon,
                "restore the default tolerance for intersection computations");
        cls.def("allowEquiDimensionalIntersections", &Constraints::allowEquiDimensionalIntersections,
                "allows/prohibits entities to intersect in geometries of the same dimension");
    }

    template<class Geo1, class Geo2, class Constraints>
    void registerBinaryEvaluator(py::class_<Constraints>& cls)
    {
        cls.def("evaluate",
                &Constraints::template evaluateBinary<Geo1, Geo2>,
                "evaluate the constraints between two entites");
    }

    template<class ctype, class Constraints>
    void registerEvaluators(py::class_<Constraints>& cls)
    {
        // types for which this ought to be able to evaluate
        using Disk = Frackit::Disk<ctype>;
        using Quad_3 = Frackit::Quadrilateral<ctype, 3>;
        using Face = OCCUtilities::FaceWrapper;

        registerBinaryEvaluator<Disk, Disk>(cls);
        registerBinaryEvaluator<Disk, Quad_3>(cls);
        registerBinaryEvaluator<Disk, Face>(cls);

        registerBinaryEvaluator<Quad_3, Quad_3>(cls);
        registerBinaryEvaluator<Quad_3, Disk>(cls);
        registerBinaryEvaluator<Quad_3, Face>(cls);

        registerBinaryEvaluator<Face, Face>(cls);
        registerBinaryEvaluator<Face, Disk>(cls);
        registerBinaryEvaluator<Face, Quad_3>(cls);
    }

} // end namespace Detail

template<class ctype>
void registerConstraints(py::module& module)
{
    using Constraints = Detail::EntityNetworkConstraintsWrapper<ctype>;
    py::class_<Constraints> cls(module, "_EntityNetworkConstraints");
    cls.def(py::init());

    Detail::registerSetters(cls);
    Detail::registerEvaluators<ctype>(cls);
}

} // end namespace Frackit::Python

#endif
