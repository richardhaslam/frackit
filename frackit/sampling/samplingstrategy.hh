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
 * \brief Interface for sampling strategies, which can be used in conjunction
 *        with, for instance, instances of 'MultiGeometrySampler'. A strategy
 *        allows for defining a number of ids, and it provides the function
 *        getNextId(), which returns an id following a specific strategy. An
 *        exemplary implementation can be found in the class
 *        'SequentialSamplingStrategy', which returns the defined ids successively.
 */
#ifndef FRACKIT_SAMPLING_STRATEGY_HH
#define FRACKIT_SAMPLING_STRATEGY_HH

#include <algorithm>
#include <vector>

#include <frackit/common/id.hh>

namespace Frackit {

/*!
 * \brief Interface for sampling strategies.
 */
class SamplingStrategy
{

public:

    /*!
     * \brief Returns the next id following the srategy.
     */
    virtual Id getNextId() = 0;

    /*!
     * \brief Adds a new id to be considered.
     */
    void addId(const Id& id)
    {
        if (!std::count(idList_.begin(), idList_.end(), id))
            idList_.push_back(id);
    }

    /*!
     * \brief Returns true if an id has been added.
     */
    bool hasId(const Id& id) const
    { return std::count(idList_.begin(), idList_.end(), id); }

protected:
    std::vector<Id> idList_;
};

} // end namespace Frackit

#endif // FRACKIT_SAMPLING_STRATEGY_HH
