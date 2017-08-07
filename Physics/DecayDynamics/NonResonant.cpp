/*
 * NonResonant.cpp
 *
 *  Created on: Jan 13, 2015
 *      Author: weidenka
 */

#include "Physics/DecayDynamics/NonResonant.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayDynamics {

std::shared_ptr<FunctionTree>
NonResonant::GetTree(const ParameterList &sample, int pos, std::string suffix) {

  int sampleSize = sample.GetMultiDouble(0)->GetNValues();
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

} /* namespace DecayDynamics */
} /* namespace Physics */
} /* namespace ComPWA */
