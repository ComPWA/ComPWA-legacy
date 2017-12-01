// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Physics/DecayDynamics/NonResonant.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayDynamics {

std::shared_ptr<FunctionTree>
NonResonant::GetTree(const ParameterList &sample, int pos, std::string suffix) {

  int sampleSize = sample.GetMultiDouble(0)->numValues();
  std::shared_ptr<MultiComplex> unitVec(
      new MultiComplex("unit", std::vector<std::complex<double>>(
                                   sampleSize, std::complex<double>(1, 0))));

  std::shared_ptr<FunctionTree> newTree(new FunctionTree());
  newTree->createHead("NonResonant" + suffix,
                      std::shared_ptr<Strategy>(new MultAll(ParType::MCOMPLEX)),
                      sampleSize);
  newTree->createLeaf("unit", unitVec, "NonResonant" + suffix); // nonReso
  return newTree;
}

} // ns::DecayDynamics
} // ns::Physics
} // ns::ComPWA
