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
 * \brief Class representing a network of entities,
 *        consisting of fragments of primary entities.
 */
#ifndef FRACKIT_ENTITY_NETWORK_HH
#define FRACKIT_ENTITY_NETWORK_HH

#include <TopTools_DataMapOfShapeInteger.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS_Shape.hxx>

#include "entitynetworkinterface.hh"

namespace Frackit {

/*!
 * \relates EntityNetworkBuilder
 * \brief Class representing a network of entities,
 *        consisting of fragments of primary entities.
 *        In addition to the fragments, it contains a
 *        map which maps an index, the index of the primary
 *        entity, to each fragment.
 * \note Use the class EntityNetworkBuilder for construction.
 */
class EntityNetwork
: public EntityNetworkInterface
{

public:
    /*!
     * \brief Constructor.
     * \param entityDim Dimension of the network entities
     * \param fragmentList List containing all entity fragments
     * \param fragmentIndexMap Map containing the fragments
     *                         of a network, where each fragment
     *                         is mapped to the index of the primary
     *                         entity of the network (before fragmentation).
     */
    EntityNetwork(int entityDim,
                  TopTools_ListOfShape&& fragmentList,
                  TopTools_DataMapOfShapeInteger&& fragmentIndexMap)
    : EntityNetworkInterface(entityDim)
    , entityFragments_(std::move(fragmentList))
    , entityFragmentIndexMap_(std::move(fragmentIndexMap))
    {}

    /*!
     * \brief Returns the fragments of the network entities.
     */
    const TopTools_ListOfShape& entityFragments() const override
    { return entityFragments_; }

    /*!
     * \brief Returns the map which maps each fragment to its primary entity index.
     */
    const TopTools_DataMapOfShapeInteger& entityFragmentsIndexMap() const override
    { return entityFragmentIndexMap_; }

private:
    TopTools_ListOfShape entityFragments_;
    TopTools_DataMapOfShapeInteger entityFragmentIndexMap_;
};

} // end namespace Frackit

#endif // FRACKIT_ENTITY_NETWORK_HH
