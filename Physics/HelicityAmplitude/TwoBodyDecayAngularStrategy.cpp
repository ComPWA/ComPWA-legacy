/*
 * TwoBodyDecayAngularStrategy.cpp
 *
 *  Created on: Aug 9, 2016
 *      Author: steve
 */

#include "Core/DataPointStorage.hpp"
#include "Physics/HelicityAmplitude/TwoBodyDecayAngularStrategy.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

TwoBodyDecayAngularStrategy::TwoBodyDecayAngularStrategy(
    std::shared_ptr<TwoBodyDecayAmplitude> tbd_amp) :
    Strategy(ParType::MCOMPLEX), tbd_amp_(tbd_amp) {
}

TwoBodyDecayAngularStrategy::~TwoBodyDecayAngularStrategy() {
  // TODO Auto-generated destructor stub
}

//! Pure Virtual interface for streaming info about the strategy
const std::string TwoBodyDecayAngularStrategy::to_str() const {
  return "WignerD";
}

//! Pure Virtual interface for executing a strategy
bool TwoBodyDecayAngularStrategy::execute(ParameterList& paras,
    std::shared_ptr<AbsParameter>& out) {

  if (ParType::MCOMPLEX) {
    std::vector<std::complex<double> > values;
    values.reserve(DataPointStorage::Instance().getNumberOfEvents());
    for (unsigned int i = 0;
        i < DataPointStorage::Instance().getNumberOfEvents(); ++i) {
      std::complex<double> value(0.0);
      auto const& evaluation_list =
          paras.GetMultiUnsignedIntegers()[0]->GetValues();
      //should be just one entry in that list
      for (unsigned int j = 0; j < evaluation_list.size(); ++j) {
        value += tbd_amp_->evaluate(i, evaluation_list[j]);
      }
      values.push_back(value);
    }
    out = std::shared_ptr<MultiComplex>(new MultiComplex(out->GetName(), values));
  }
  return true;
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
