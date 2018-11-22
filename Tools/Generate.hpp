// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Some useful function for Monte-Carlo generation.
///

#ifndef COMPWA_TOOLS_GENERATE_HPP_
#define COMPWA_TOOLS_GENERATE_HPP_

#include "Core/AmpIntensity.hpp"
#include "Core/Generator.hpp"
#include "Core/Kinematics.hpp"
#include "Core/ProgressBar.hpp"
#include "Data/Data.hpp"
#include "Tools/Integration.hpp"

namespace ComPWA {
namespace Tools {

std::shared_ptr<ComPWA::Data::Data>
generate(unsigned int NumberOfEvents,
         std::shared_ptr<ComPWA::Kinematics> Kinematics,
         std::shared_ptr<ComPWA::Generator> Generator,
         std::shared_ptr<ComPWA::AmpIntensity> Intensity);

std::shared_ptr<ComPWA::Data::Data>
generate(unsigned int NumberOfEvents,
         std::shared_ptr<ComPWA::Kinematics> Kinematics,
         std::shared_ptr<ComPWA::Generator> Generator,
         std::shared_ptr<ComPWA::AmpIntensity> Intensity,
         std::shared_ptr<ComPWA::Data::Data> phsp,
         std::shared_ptr<ComPWA::Data::Data> phspTrue =
             std::shared_ptr<ComPWA::Data::Data>());

inline std::shared_ptr<ComPWA::Data::Data>
generatePhsp(unsigned int nEvents, std::shared_ptr<ComPWA::Generator> gen) {
  std::shared_ptr<ComPWA::Data::Data> sample(new ComPWA::Data::Data);

  LOG(INFO) << "Generating phase-space MC: [" << nEvents << " events] ";

  ComPWA::ProgressBar bar(nEvents);
  for (unsigned int i = 0; i < nEvents; ++i) {
    ComPWA::Event tmp = gen->generate();
    double ampRnd = gen->uniform(0, 1);
    if (ampRnd > tmp.weight()) {
      --i;
      continue;
    }

    // Reset weights: weights are taken into account by hit&miss. The
    // resulting sample is therefore unweighted
    tmp.setWeight(1.);

    tmp.setEfficiency(1.);
    sample->add(tmp);
    bar.next();
  }
  return sample;
}

std::shared_ptr<ComPWA::Data::Data>
generateImportanceSampledPhsp(unsigned int NumberOfEvents,
                              std::shared_ptr<ComPWA::Kinematics> Kinematics,
                              std::shared_ptr<ComPWA::Generator> Generator,
                              std::shared_ptr<ComPWA::AmpIntensity> Intensity);

} // namespace Tools
} // namespace ComPWA
#endif
