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
 * \brief Sampling strategy implementation which goes
 *        over the defined ids successively.
 */
#ifndef FRACKIT_SEQUENTIAL_SAMPLING_STRATEGY_HH
#define FRACKIT_SEQUENTIAL_SAMPLING_STRATEGY_HH

#include <stdexcept>

#include "samplingstrategy.hh"

namespace Frackit {

/*!
 * \brief Sampling strategy implementation which goes
 *        over the defined ids successively.
 */
class SequentialSamplingStrategy
: public SamplingStrategy
{

public:

    /*!
     * \brief Returns the next id.
     */
    Id getNextId() override
    {
        if (this->idList_.empty())
            throw std::runtime_error("Id list is empty!");

        // The id about to returned
        const auto curId = this->idList_[curIdx_];

        // Update internal index for next call
        if (curIdx_ != this->idList_.size() - 1) curIdx_++;
        else curIdx_ = 0;

        return curId;
    }

private:
    std::size_t curIdx_ = 0;
};

} // end namespace Frackit

#endif // FRACKIT_SEQUENTIAL_SAMPLING_STRATEGY_HH
