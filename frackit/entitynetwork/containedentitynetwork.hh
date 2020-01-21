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
  * \brief Class representing a network of entities, contained
  *        in (possibly multiple) sub-domains. Sub-networks might
  *        be defined on each sub-domain.
  */
#ifndef FRACKIT_CONTAINED_ENTITY_NETWORK_HH
#define FRACKIT_CONTAINED_ENTITY_NETWORK_HH

#include <vector>
#include <unordered_map>

#include <TopTools_DataMapOfShapeInteger.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS_Shape.hxx>

#include <frackit/common/id.hh>
#include "containedentitynetworkinterface.hh"

namespace Frackit {

/*!
 * \relates ContainedEntityNetworkBuilder
 * \brief Class representing a network of entities, contained
 *        in (possibly multiple) sub-domains. Sub-networks might
 *        be defined on each sub-domain.
 * \note Use the class ContainedEntityNetworkBuilder for construction.
 */
class ContainedEntityNetwork
: public ContainedEntityNetworkInterface
{

public:
    /*!
     * \brief Constructor.
     * \param entityDim Dimension of the network entities
     * \param fragmentList List containing all entity fragments
     * \param fragmentIndexMap Map containing the fragments
     *                         of a network, where each fragment
     *                         is mapped to the index of the primary
     *                         entity of the network.
     */
    ContainedEntityNetwork(int entityDim,
                           int domainDim,
                           std::unordered_map<std::size_t, TopTools_ListOfShape>&& sdFragments,
                           std::unordered_map<std::size_t, TopTools_ListOfShape>&& entityFragments,
                           std::unordered_map<std::size_t, TopTools_DataMapOfShapeInteger>&& entityFragmentMaps)
    : ContainedEntityNetworkInterface(entityDim, domainDim)
    , subDomainFragments_(std::move(sdFragments))
    , subDomainEntityFragments_(std::move(entityFragments))
    , subDomainEntityFragmentIndexMap_(std::move(entityFragmentMaps))
    {
        subDomainIds_.reserve(subDomainFragments_.size());
        for (const auto& sdDataPair : subDomainFragments_)
            subDomainIds_.emplace_back(sdDataPair.first);
    }

    /*!
     * \brief Returns the ids of defined the sub-domains
     */
    const std::vector<Id>& subDomainIds() const override
    { return subDomainIds_; }

    /*!
     * \brief Returns the fragments of a sub-domain
     * \param subDomainIdx The index of the sub-domain
     */
    const TopTools_ListOfShape& subDomainFragments(Id subDomainId) const override
    { return subDomainFragments_.at(subDomainId.get()); }

    /*!
     * \brief Returns the entity fragments of the network defined for a sub-domain
     * \param subDomainIdx The index of the sub-domain
     */
    const TopTools_ListOfShape& subDomainEntityFragments(Id subDomainId) const override
    { return subDomainEntityFragments_.at(subDomainId.get()); }

    /*!
     * \brief Returns the map which maps each fragment the network of a sub-domain to its primary entity index.
     * \param subDomainIdx The index of the sub-domain
     */
    const TopTools_DataMapOfShapeInteger& subDomainEntityFragmentsIndexMap(Id subDomainId) const override
    { return subDomainEntityFragmentIndexMap_.at(subDomainId.get()); }

private:
    std::unordered_map<std::size_t, TopTools_ListOfShape> subDomainFragments_;
    std::unordered_map<std::size_t, TopTools_ListOfShape> subDomainEntityFragments_;
    std::unordered_map<std::size_t, TopTools_DataMapOfShapeInteger> subDomainEntityFragmentIndexMap_;
    std::vector<Id> subDomainIds_;
};

} // end namespace Frackit

#endif // FRACKIT_CONTAINED_ENTITY_NETWORK_HH
