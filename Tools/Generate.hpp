// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Some useful functions for Monte-Carlo event generation.
///

#ifndef COMPWA_TOOLS_GENERATE_HPP_
#define COMPWA_TOOLS_GENERATE_HPP_

#include <memory>

namespace ComPWA {

class Kinematics;
class Generator;
class Intensity;

namespace Data {
class DataSet;
}

namespace Tools {

std::shared_ptr<ComPWA::Data::DataSet>
generate(unsigned int NumberOfEvents,
         std::shared_ptr<ComPWA::Kinematics> Kinematics,
         std::shared_ptr<ComPWA::Generator> Generator,
         std::shared_ptr<ComPWA::Intensity> Intensity);

std::shared_ptr<ComPWA::Data::DataSet>
generate(unsigned int NumberOfEvents,
         std::shared_ptr<ComPWA::Kinematics> Kinematics,
         std::shared_ptr<ComPWA::Generator> Generator,
         std::shared_ptr<ComPWA::Intensity> Intensity,
         std::shared_ptr<ComPWA::Data::DataSet> phsp,
         std::shared_ptr<ComPWA::Data::DataSet> phspTrue = {});

std::shared_ptr<ComPWA::Data::DataSet>
generatePhsp(unsigned int nEvents, std::shared_ptr<ComPWA::Generator> gen);

std::shared_ptr<ComPWA::Data::DataSet>
generateImportanceSampledPhsp(unsigned int NumberOfEvents,
                              std::shared_ptr<ComPWA::Kinematics> Kinematics,
                              std::shared_ptr<ComPWA::Generator> Generator,
                              std::shared_ptr<ComPWA::Intensity> Intensity);

} // namespace Tools
} // namespace ComPWA
#endif
