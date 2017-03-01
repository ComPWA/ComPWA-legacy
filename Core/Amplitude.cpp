/*
 * Amplitude.cpp
 *
 *  Created on: Mar 16, 2016
 *      Author: weidenka
 */

#include "Core/Amplitude.hpp"

namespace ComPWA {

void Amplitude::UpdateParameters(ParameterList &par) {
  //TODO:

  return;
}

void Amplitude::GetParameters(ParameterList &outPar) const {
  outPar.AddParameter(_magnitude);
  outPar.AddParameter(_phase);
}

} /* namespace ComPWA */
