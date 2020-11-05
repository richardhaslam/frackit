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
#ifndef FRACKIT_PYTHON_INTERSECT_HH
#define FRACKIT_PYTHON_INTERSECT_HH

#include <vector>
#include <variant>
#include <type_traits>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>

#include <frackit/geometry/segment.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/polygon.hh>
#include <frackit/geometry/cylindersurface.hh>

#include <frackit/python/occutilities/brepwrapper.hh>
#include <frackit/intersection/intersect.hh>
#include <frackit/intersection/intersectiontraits.hh>

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

//! for compatibility with BRepWrapper classes
template<class Geo> struct UnwrapperHelper { using type = Geo; };
template<class Geo> struct UnwrapperHelper<OCCUtilities::BRepWrapper<Geo>> { using type = Geo; };
template<class Geo> using UnwrappedType = typename UnwrapperHelper<Geo>::type;

//! the converted type of an individual geometry type
template<class Geo>
using ConvertedType = std::conditional_t< std::is_convertible_v<Geo, TopoDS_Shape>,
                                          OCCUtilities::BRepWrapper<Geo>,
                                          Geo >;

template<class T> struct IsVector : public std::false_type {};
template<class T> struct IsVector<std::vector<T>> : public std::true_type {};

template<class T> struct IntersectionVariant { using type = T; };
template<class T> struct IntersectionVariant<std::vector<T>> { using type = T; };

template<class T> struct ConvertedVariant;
template<class... Geoms>
struct ConvertedVariant<std::variant<Geoms...>>
{ using type = std::variant< ConvertedType<Geoms>... >; };

//! Convert an intersection result into new one using ShapeWrappers
template<class T>
using ConvertVariant = typename ConvertedVariant<T>::type;

//! Uses the converter to define the new intersection result
template<class T>
using ConvertIntersection = std::conditional_t< IsVector<T>::value,
                                                std::vector<ConvertVariant<typename IntersectionVariant<T>::type>>,
                                                ConvertVariant<typename IntersectionVariant<T>::type> >;

//! Alias for the converted intersection result type
template<class Intersection>
using ConvertedIntersection = ConvertIntersection<Intersection>;

//! Alias for the intersection result type between two geometries
template<class Geo1, class Geo2>
using IntersectionResult = ConvertedIntersection<Frackit::Intersection<UnwrappedType<Geo1>,
                                                                       UnwrappedType<Geo2>>>;

//! convert shape into wrapper
template<class Geo, std::enable_if_t<std::is_convertible_v<Geo, TopoDS_Shape>, int> = 0>
OCCUtilities::BRepWrapper<Geo> convertShape(const Geo& geo) { return {geo}; }

template<class Geo, std::enable_if_t<!std::is_convertible_v<Geo, TopoDS_Shape>, int> = 0>
Geo convertShape(const Geo& geo) { return geo; }

//! parse the result of an intersection into the converted type
template<class Intersection, std::enable_if_t<!IsVector<Intersection>::value, int> = 0>
ConvertedIntersection<Intersection> convertIntersectionResult(const Intersection& is)
{
    ConvertedIntersection<Intersection> result;
    std::visit([&result] (const auto& geo) { result = convertShape(geo); }, is);
    return result;
 }

//! parse the result of an intersection (multiple geometries) into the converted type
template<class Intersection, std::enable_if_t<IsVector<Intersection>::value, int> = 0>
ConvertedIntersection<Intersection> convertIntersectionResult(const Intersection& is)
{
    ConvertedIntersection<Intersection> result;
    for (const auto& geo : is)
        result.push_back(convertIntersectionResult(geo));
    return result;
}

//! convert intersection result such that all TopoDS_Shape classes
//! are converted into the ShapeWrappers registered in the Python bindings
template<class Geo1, class Geo2>
IntersectionResult<Geo1, Geo2> intersectAndConvert(const Geo1& geo1, const Geo2& geo2)
{ return convertIntersectionResult(Frackit::intersect(geo1, geo2)); }

//! convert intersection result such that all TopoDS_Shape classes
//! are converted into the ShapeWrappers registered in the Python bindings
template<class Geo1, class Geo2, class ctype>
IntersectionResult<Geo1, Geo2> intersectAndConvert(const Geo1& geo1, const Geo2& geo2, ctype eps)
{ return convertIntersectionResult(Frackit::intersect(geo1, geo2, eps)); }

//! convert shape wrapper function arguments into shape representations and forward
template<class Geo1, class Geo2>
IntersectionResult<Geo1, Geo2> intersect(const Geo1& geo1, const Geo2& geo2)
{
    using namespace OCCUtilities;
    return intersectAndConvert(getUnwrappedShape(geo1), getUnwrappedShape(geo2));
}

//! convert shape wrapper function arguments into shape representations and forward
template<class Geo1, class Geo2, class ctype>
IntersectionResult<Geo1, Geo2> intersect(const Geo1& geo1, const Geo2& geo2, ctype eps)
{
    using namespace OCCUtilities;
    return intersectAndConvert(getUnwrappedShape(geo1), getUnwrappedShape(geo2), eps);
}


template<class Geo1, class Geo2, class ctype>
void registerIntersectionFunction(py::module& module,
                                  const std::string& name1,
                                  const std::string& name2)
{
    using namespace py::literals;

    const std::string doc = "Returns the intersection between the geometries " + name1 + "/" + name2;
    const std::string docDef = doc + " (default geometric tolerance)";
    const std::string docGiv = doc + " (given geometric tolerance)";
    module.def("intersect",
               py::overload_cast<const Geo1&, const Geo2&>(&Detail::intersect<Geo1, Geo2>),
               docDef.c_str(), "geo1"_a, "geo2"_a);
    module.def("intersect",
               py::overload_cast<const Geo1&, const Geo2&, ctype>(&Detail::intersect<Geo1, Geo2, ctype>),
               docGiv.c_str(), "geo1"_a, "geo2"_a, "eps"_a);
}

} // end namespace Detail

template<class ctype>
void registerIntersectionFunctions(py::module& module)
{
    using Segment = Frackit::Segment<ctype, 3>;
    using Quad = Frackit::Quadrilateral<ctype, 3>;
    using Polygon = Frackit::Polygon<ctype, 3>;
    using Disk = Frackit::Disk<ctype>;
    using CylSurf = Frackit::CylinderSurface<ctype>;
    using FaceWrapper = Frackit::Python::OCCUtilities::FaceWrapper;

    Detail::registerIntersectionFunction<Segment, Segment, ctype>(module, "Segment", "Segment");
    Detail::registerIntersectionFunction<Quad, Quad, ctype>(module, "Quadrilateral", "Quadrilateral");
    Detail::registerIntersectionFunction<Polygon, Polygon, ctype>(module, "Polygon", "Polygon");
    Detail::registerIntersectionFunction<Disk, Disk, ctype>(module, "Disk", "Disk");
    Detail::registerIntersectionFunction<FaceWrapper, FaceWrapper, ctype>(module, "TopoDS_Face_Wrapper", "TopoDS_Face_Wrapper");

    Detail::registerIntersectionFunction<Quad, Disk, ctype>(module, "Quadrilateral", "Disk");
    Detail::registerIntersectionFunction<Disk, Quad, ctype>(module, "Disk", "Quadrilateral");

    Detail::registerIntersectionFunction<Quad, Disk, ctype>(module, "Polygon", "Disk");
    Detail::registerIntersectionFunction<Disk, Quad, ctype>(module, "Disk", "Polygon");

    Detail::registerIntersectionFunction<FaceWrapper, Disk, ctype>(module, "TopoDS_Face_Wrapper", "Disk");
    Detail::registerIntersectionFunction<Disk, FaceWrapper, ctype>(module, "Disk", "TopoDS_Face_Wrapper");

    Detail::registerIntersectionFunction<FaceWrapper, Quad, ctype>(module, "TopoDS_Face_Wrapper", "Quadrilateral");
    Detail::registerIntersectionFunction<Quad, FaceWrapper, ctype>(module, "Quadrilateral", "TopoDS_Face_Wrapper");

    Detail::registerIntersectionFunction<FaceWrapper, Polygon, ctype>(module, "TopoDS_Face_Wrapper", "Polygon");
    Detail::registerIntersectionFunction<Polygon, FaceWrapper, ctype>(module, "Polygon", "TopoDS_Face_Wrapper");

    Detail::registerIntersectionFunction<CylSurf, Disk, ctype>(module, "CylinderSurface", "Disk");
    Detail::registerIntersectionFunction<Disk, CylSurf, ctype>(module, "Disk", "CylinderSurface");

    Detail::registerIntersectionFunction<CylSurf, Quad, ctype>(module, "CylinderSurface", "Quadrilateral");
    Detail::registerIntersectionFunction<Quad, CylSurf, ctype>(module, "Quadrilateral", "CylinderSurface");

    Detail::registerIntersectionFunction<CylSurf, Polygon, ctype>(module, "CylinderSurface", "Polygon");
    Detail::registerIntersectionFunction<Polygon, CylSurf, ctype>(module, "Polygon", "CylinderSurface");

    Detail::registerIntersectionFunction<CylSurf, FaceWrapper, ctype>(module, "CylinderSurface", "TopoDS_Face_Wrapper");
    Detail::registerIntersectionFunction<FaceWrapper, CylSurf, ctype>(module, "TopoDS_Face_Wrapper", "CylinderSurface");
}

} // end namespace Frackit::Python

#endif
