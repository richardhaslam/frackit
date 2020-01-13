// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
/*!
 * \file
 * \brief \todo TODO doc me.
 */
#ifndef FRACKIT_GP_UTILITIES_HH
#define FRACKIT_GP_UTILITIES_HH

// objects from geometric processors package
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <gp_Vec.hxx>

// internal geometry classes
#include <frackit/geometry/point.hh>
#include <frackit/geometry/direction.hh>
#include <frackit/geometry/vector.hh>

namespace Frackit {
namespace OCCUtilities {

    //! converts an n-d point into a 3d point
    template<class ctype, int dim>
    Point<ctype, 3> convertTo3d(const Point<ctype, dim>& p)
    {
        static_assert(dim <= 3 && dim != 0, "Only 0 < dim <= 3 supported");
        if (dim == 1) return Point<ctype, 3>(p.x(), 0.0, 0.0);
        else if (dim == 2) return Point<ctype, 3>(p.x(), p.y(), 0.0);
        else if (dim == 3) return p;
    }

    //! converts a point to an object from the geometric processors package
    template<class ctype, int dim>
    gp_Pnt point(const Point<ctype, dim>& p)
    {
        static_assert(dim <= 3 && dim != 0, "Only 0 < dim <= 3 supported");
        if (dim == 1) return gp_Pnt(p.x(), 0.0, 0.0);
        else if (dim == 2) return gp_Pnt(p.x(), p.y(), 0.0);
        else if (dim == 3) return gp_Pnt(p.x(), p.y(), p.z());
    }

    //! casts a point from the geometric processors package
    Point<double, 3> point(const gp_Pnt& p)
    { return {p.X(), p.Y(), p.Z()}; }

    //! converts an internal direction into gp direction object
    template<class ctype, int dim>
    gp_Dir direction(const Direction<ctype, dim>& dir)
    {
        static_assert(dim <= 3 && dim != 0, "Only 0 < dim <= 3 supported");
        if (dim == 1) return gp_Dir(dir.x(), 0.0, 0.0);
        else if (dim == 2) return gp_Dir(dir.x(), dir.y(), 0.0);
        else if (dim == 3) return gp_Dir(dir.x(), dir.y(), dir.z());
    }

    //! converts a gp direction object into an internal one
    template<class ctype = double, int dim = 3>
    Direction<ctype, dim> direction(const gp_Dir& dir)
    {
        static_assert(dim == 3, "Currently only dim == 3 supported");
        using Vector = Frackit::Vector<ctype, dim>;
        return Direction<ctype, dim>(Vector(dir.X(), dir.Y(), dir.Z()));
    }

    //! converts a gp vector object into an internal one
    template<class ctype = double, int dim = 3>
    Vector<ctype, dim> vector(const gp_Vec& v)
    {
        static_assert(dim == 3, "Currently only dim == 3 supported");
        return Vector<ctype, dim>(v.X(), v.Y(), v.Z());
    }

} // end namespace OCCUtilities
} // end namespace Frackit

#endif // FRACKIT_GP_UTILITIES_HH
