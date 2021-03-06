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

#include <iomanip>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <utility>
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
        else if (it->second > targetCount)
            std::cout << "Warning: given target count is below current count" << std::endl;
    }

    /*!
     * \brief Reset everything.
     */
    void reset()
    {
        headerPrinted_ = false;
        rejectedCount_ = 0;
        count_.clear();
        targetCount_.clear();
    }

    /*!
     * \brief Reset all counters.
     */
    void resetCounters()
    {
        for (auto& count : count_) count.second = 0;
    }

    /*!
     * \brief Reset counter for the given id.
     */
    void resetCounter(const Id& id)
    {
        count_[id.get()] = 0;
    }

    /*!
     * \brief Returns true when the target
     *        sample count is reached.
     */
    bool finished()
    {
        for (const auto& pair : count_)
            if (pair.second < targetCount_.at(pair.first))
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
        return it->second >= targetCount_.at(it->first);
    }

    /*!
     * \brief Increase counter for a specific id.
     */
    void increaseCounter(const Id& id)
    {
        count_[id.get()]++;
        if (count_[id.get()] > targetCount_[id.get()])
            std::cout << "Warning: target count for id " << id.get() << " was surpassed" << std::endl;
    }

    /*!
     * \brief Increase counter of rejected samples.
     */
    void increaseRejectedCounter()
    {
        rejectedCount_++;
    }

    /*!
     * \brief Returns the entity count for the given id.
     */
    std::size_t getCount(const Id& id)
    { return count_[id.get()]; }

    /*!
     * \brief Returns the overall entity count.
     */
    std::size_t getCount()
    {
        return std::accumulate(count_.begin(),
                               count_.end(),
                               0,
                               [] (const auto& curCount, const auto& idCountPair)
                               { return curCount + idCountPair.second; });
    }

    /*!
     * \brief Print current status to terminal.
     */
    void print(bool forceHeaderPrint = false)
    {
        if (!headerPrinted_ || forceHeaderPrint)
        {
            std::cout << "#####################################################################################\n"
                      << "# Accepted count   |   Rejected count   |   Acceptance Ratio [%]   |   Progress [%] #\n"
                      << "#-----------------------------------------------------------------------------------#\n";
            headerPrinted_ = true;
        }

        std::size_t curCount = 0;
        std::size_t curTargetCount = 0;
        for (const auto& pair : count_) curCount += pair.second;
        for (const auto& pair : targetCount_) curTargetCount += pair.second;

        const auto ratio = 100.0*double(double(curCount)/double(curCount+rejectedCount_));
        const auto progress = 100.0*double(curCount)/double(curTargetCount);

        std::cout << std::setprecision(2) << std::fixed;
        const auto ratioNumChars = std::to_string(int(ratio)).size() + 3;
        const auto progressNumChars = std::to_string(int(progress)).size() + 3;

        const auto countString = std::to_string(curCount);
        const auto rejectedCountString = std::to_string(rejectedCount_);

        using std::max;
        const std::size_t zero = 0;
        const std::size_t paddingCount    = max(zero, 13 - countString.size());
        const std::size_t paddingRejected = max(zero, 13 - rejectedCountString.size());
        const std::size_t paddingRatio    = max(zero, 16 - ratioNumChars);
        const std::size_t paddingProgess  = max(zero, 8  - progressNumChars);

        std::cout << "  "      << std::string( paddingCount, ' ')    << countString + ' '
                  << "   |   " << std::string( paddingRejected, ' ') << rejectedCountString + ' '
                  << "   |   " << std::string( paddingRatio, ' ')    << ratio << std::string(4, ' ')
                  << "   |   " << std::string( paddingProgess, ' ')  << progress << std::endl;
    }

private:
    bool headerPrinted_ = false;
    std::size_t rejectedCount_ = 0;
    std::unordered_map<std::size_t, std::size_t> count_;
    std::unordered_map<std::size_t, std::size_t> targetCount_;
};

} // end namespace Frackit

#endif // FRACKIT_SAMPLING_STATUS_HH
