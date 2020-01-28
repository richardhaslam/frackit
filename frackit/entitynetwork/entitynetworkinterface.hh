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
/*!
 * \file
 * \ingroup EntityNetwork
 * \brief Interface defining networks of entities.
 */
#ifndef FRACKIT_ENTITY_NETWORK_INTERFACE_HH
#define FRACKIT_ENTITY_NETWORK_INTERFACE_HH

#include <TopTools_DataMapOfShapeInteger.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS_Shape.hxx>

namespace Frackit {

/*!
 * \ingroup EntityNetwork
 * \brief Interface defining networks of entities.
 */
class EntityNetworkInterface
{
public:
    /*!
     * \brief Constructor.
     * \param entityDim Dimension of the network entities
     */
    EntityNetworkInterface(int entityDim)
    : entityDimension_(entityDim)
    {}

    /*!
     * \brief Returns the dimension of the network entities.
     */
    int entityDimension() const
    { return entityDimension_; }

    /*!
     * \brief Returns the fragments of the network entities.
     */
    virtual const TopTools_ListOfShape& entityFragments() const = 0;

    /*!
     * \brief Returns the map which maps each fragment to its primary entity index.
     */
    virtual const TopTools_DataMapOfShapeInteger& entityFragmentsIndexMap() const = 0;

private:
    int entityDimension_;
};

} // end namespace Frackit

#endif // FRACKIT_ENTITY_NETWORK_INTERFACE_HH
