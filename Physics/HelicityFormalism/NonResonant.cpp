/*
 * NonResonant.cpp
 *
 *  Created on: Jan 13, 2015
 *      Author: weidenka
 */

#include "Physics/HelicityFormalism/NonResonant.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

std::shared_ptr<FunctionTree> NonResonant::GetTree(ParameterList &sample,
                                                   ParameterList &toySample,
                                                   std::string suffix) {
  int sampleSize = sample.GetMultiDouble(0)->GetNValues();
  double phspVol = Kinematics::Instance()->GetPhspVolume();

  std::shared_ptr<FunctionTree> newTree(new FunctionTree());
  newTree->createHead("DynamicalFunction",
                      std::shared_ptr<Strategy>(new MultAll(ParType::MCOMPLEX)),
                      sampleSize);
  std::shared_ptr<MultiComplex> unitVec(
      new MultiComplex("unit", std::vector<std::complex<double>>(
                                   sampleSize, std::complex<double>(1, 0))));
  newTree->createLeaf("NonResonant" + _name, unitVec,
                      "DynamicalFunction"); // nonReso
  newTree->createLeaf("Normalization" + _name, 1 / std::sqrt(phspVol),
                      "DynamicalFunction");
  return newTree;
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
