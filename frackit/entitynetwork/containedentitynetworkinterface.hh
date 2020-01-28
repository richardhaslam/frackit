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
 * \brief Interface defining networks of entities contained
 *        in (possibly multiple) sub-domains.
 */
#ifndef FRACKIT_CONTAINED_ENTITY_NETWORK_INTERFACE_HH
#define FRACKIT_CONTAINED_ENTITY_NETWORK_INTERFACE_HH

#include <vector>

#include <TopTools_DataMapOfShapeInteger.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS_Shape.hxx>

#include <frackit/common/id.hh>

namespace Frackit {

/*!
 * \ingroup EntityNetwork
 * \brief Interface defining networks of entities contained
 *        in (possibly multiple) sub-domains.
 */
class ContainedEntityNetworkInterface
{

public:

    /*!
     * \brief Constructor.
     * \param entityDim Dimension of the network entities
     * \param domainDim Dimension of the domain in
     *                  which entities are embedded
     */
    ContainedEntityNetworkInterface(int entityDim, int domainDim)
    : entityDimension_(entityDim)
    , domainDimension_(domainDim)
    {}

    /*!
     * \brief Returns the dimension of the network entities.
     */
    int entityDimension() const
    { return entityDimension_; }

    /*!
     * \brief Returns the dimension of the (sub-)domains.
     */
    int domainDimension() const
    { return domainDimension_; }

    /*!
     * \brief Returns the ids of defined the sub-domains
     */
    virtual const std::vector<Id>& subDomainIds() const = 0;

    /*!
     * \brief Returns the fragments of a sub-domain
     * \param subDomainId The id of the sub-domain
     */
    virtual const TopTools_ListOfShape& subDomainFragments(Id subDomainId) const = 0;

    /*!
     * \brief Returns the entity fragments of the network defined for a sub-domain
     * \param subDomainId The id of the sub-domain
     */
    virtual const TopTools_ListOfShape& subDomainEntityFragments(Id subDomainId) const = 0;

    /*!
     * \brief Returns the map which maps each fragment the network of a sub-domain to its primary entity index.
     * \param subDomainId The id of the sub-domain
     */
    virtual const TopTools_DataMapOfShapeInteger& subDomainEntityFragmentsIndexMap(Id subDomainId) const = 0;

private:
    int entityDimension_;
    int domainDimension_;
};

} // end namespace Frackit

#endif // FRACKIT_CONTAINED_ENTITY_NETWORK_INTERFACE_HH
