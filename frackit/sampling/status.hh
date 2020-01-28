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
 * \ingroup Sampling
 * \brief Class that can be used to define target counters in
 *        sampling procedures and which outputs the current status.
 */
#ifndef FRACKIT_SAMPLING_STATUS_HH
#define FRACKIT_SAMPLING_STATUS_HH

#include <unordered_map>
#include <initializer_list>

#include <frackit/common/id.hh>

namespace Frackit {

/*!
 * \ingroup Sampling
 * \brief Class that can be used to define target counters in
 *        sampling procedures and which outputs the current status.
 */
class SamplingStatus
{

public:

    /*!
     * \brief Set the target sample count for a specific id.
     */
    void setTargetCount(const Id& id, std::size_t targetCount)
    {
        targetCount_[id.get()] = targetCount;

        auto it = count_.find(id.get());
        if (it == count_.end())
            count_[id.get()] = 0;
    }

    /*!
     * \brief Returns true when the target
     *        sample count is reached.
     */
    bool finished()
    {
        for (const auto& pair : count_)
            if (pair.second != targetCount_.at(pair.first))
                return false;
        return true;
    }

    /*!
     * \brief Returns true when the target count
     *        for a specific id is reached.
     */
    bool finished(const Id& id)
    {
        auto it = count_.find(id.get());
        if (it == count_.end())
            throw std::runtime_error("Target count not set for given id");
        return it->second == targetCount_.at(it->first);
    }

    /*!
     * \brief Increase counter for a specific id.
     */
    void increaseCounter(const Id& id)
    {
        count_[id.get()]++;
    }

    /*!
     * \brief Increase counter of rejected samples.
     */
    void increaseRejectedCounter()
    {
        rejectedCount_++;
    }

    /*!
     * \brief Print current status to terminal.
     */
    void print(bool forceHeaderPrint = false)
    {
        if (!headerPrinted_ || forceHeaderPrint)
        {
            std::cout << "#################################################################################\n"
                      << "# Accepted count   |   Rejected count   |   Acceptance Ratio   |   Progress [%] #\n"
                      << "#-------------------------------------------------------------------------------#\n";
            headerPrinted_ = true;
        }

        std::size_t curCount = 0;
        std::size_t curTargetCount = 0;
        for (const auto& pair : count_) curCount += pair.second;
        for (const auto& pair : targetCount_) curTargetCount += pair.second;

        const auto ratio = double(double(curCount)/double(curCount+rejectedCount_));
        const auto progress = 100.0*double(curCount)/double(curTargetCount);

        std::cout << "  "   << std::setw(17) << std::setfill(' ')
                            << std::to_string(curCount) + std::string(9, ' ')
                  << "|   " << std::setprecision(6) << std::setw(17) << std::setfill(' ')
                            << std::to_string(rejectedCount_) + std::string(9, ' ')
                  << "|   " << std::setprecision(6) << std::setw(19) << std::setfill(' ')
                            << std::to_string( ratio ) + std::string(7, ' ');
        std::cout << "|   " << std::setprecision(2) << std::string(2, ' ') + std::to_string( progress ) << std::endl;
    }

private:
    bool headerPrinted_ = false;
    std::size_t rejectedCount_ = 0;
    std::unordered_map<std::size_t, std::size_t> count_;
    std::unordered_map<std::size_t, std::size_t> targetCount_;
};

} // end namespace Frackit

#endif // FRACKIT_SAMPLING_STATUS_HH
