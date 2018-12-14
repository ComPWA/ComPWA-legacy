// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Some useful functions for Monte-Carlo event generation.
///

#ifndef COMPWA_TOOLS_GENERATE_HPP_
#define COMPWA_TOOLS_GENERATE_HPP_

#include "Core/Generator.hpp"
#include "Core/Intensity.hpp"
#include "Core/Kinematics.hpp"
#include "Core/ProgressBar.hpp"

namespace ComPWA {
namespace Tools {

std::vector<ComPWA::Event>
generate(unsigned int NumberOfEvents,
         std::shared_ptr<ComPWA::Kinematics> Kinematics,
         std::shared_ptr<ComPWA::Generator> Generator,
         std::shared_ptr<ComPWA::Intensity> Intensity);

std::vector<ComPWA::Event> generate(
    unsigned int NumberOfEvents, std::shared_ptr<ComPWA::Kinematics> Kinematics,
    std::shared_ptr<ComPWA::Generator> Generator,
    std::shared_ptr<ComPWA::Intensity> Intensity,
    const std::vector<ComPWA::Event> &phsp,
    const std::vector<ComPWA::Event> &phspTrue = std::vector<ComPWA::Event>());

inline std::vector<ComPWA::Event>
generatePhsp(unsigned int nEvents, std::shared_ptr<ComPWA::Generator> gen) {
  std::vector<ComPWA::Event> sample;

  LOG(INFO) << "Generating phase-space MC: [" << nEvents << " events] ";

  ComPWA::ProgressBar bar(nEvents);
  for (unsigned int i = 0; i < nEvents; ++i) {
    ComPWA::Event tmp = gen->generate();
    double ampRnd = gen->uniform(0, 1);
    if (ampRnd > tmp.Weight) {
      --i;
      continue;
    }

    // Reset weights: weights are taken into account by hit&miss. The
    // resulting sample is therefore unweighted
    tmp.Weight = 1.0;

    tmp.Efficiency = 1.0;
    sample.push_back(tmp);
    bar.next();
  }
  return sample;
}

std::vector<ComPWA::Event>
generateImportanceSampledPhsp(unsigned int NumberOfEvents,
                              std::shared_ptr<ComPWA::Kinematics> Kinematics,
                              std::shared_ptr<ComPWA::Generator> Generator,
                              std::shared_ptr<ComPWA::Intensity> Intensity);

} // namespace Tools
} // namespace ComPWA
#endif
