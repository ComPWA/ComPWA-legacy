// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_DATA_DATATRANSFORMATION_HPP_
#define COMPWA_DATA_DATATRANSFORMATION_HPP_

#include <vector>

#include "Core/Event.hpp"
#include "Core/ParameterList.hpp"

namespace ComPWA {
class Kinematics;
namespace Data {

std::vector<DataPoint>
convertEventsToDataPoints(const std::vector<Event> &Events,
                             std::shared_ptr<ComPWA::Kinematics> Kinematics);

ParameterList
convertEventsToParameterList(const std::vector<Event> &Events,
                             std::shared_ptr<ComPWA::Kinematics> Kinematics);

void reduceToPhaseSpace(std::vector<Event> &Events,
                        std::shared_ptr<ComPWA::Kinematics> Kinematics);

/// Get 'horizontal' list of the kinematic variables. For each variable
/// (e.g. m23sq, m13sq ...) a MultiDouble is added to ParameterList.
/// This ParameterList is used in the FunctionTree.
ComPWA::ParameterList
convertDataPointsToParameterList(const std::vector<DataPoint> &DataPointList);

std::vector<DataPoint> convertParameterListToDataPoints(
    const ComPWA::ParameterList &DataParameterList);

/// Select random sub sample of data sets.
/// A hit&miss procedure is applied to the first sample to select a random set
/// of events. In case an event of the first sample is added to the output
/// sample, an event from the second sample is also added to the corresponding
/// output sample. This function can be used for two sample that are in sync.
/// E.g. a sample with reconstructed values and a sample with the
/// corresponding true values We expect that @param out1 and @param out2 are
/// pointers to empty samples.
///
/// \param size size of sub sample
/// \param gen generator
/// \param in1 input sample 1
/// \param in2 input sample 2
/// \param out1 output sub sample from sample 1
/// \param out2 output sub sample from sample 2
/* void rndReduceSet(unsigned int size, std::shared_ptr<Generator> gen,
                         Data *in1, Data *out1, Data *in2 = NULL,
                         Data *out2 = NULL);*/

} // namespace Data
} // namespace ComPWA

#endif
