/*
 * TwoBodyDecayAngularStrategy.hpp
 *
 *  Created on: Aug 9, 2016
 *      Author: steve
 */

#ifndef PHYSICS_HELICITYAMPLITUDE_TWOBODYDECAYANGULARSTRATEGY_HPP_
#define PHYSICS_HELICITYAMPLITUDE_TWOBODYDECAYANGULARSTRATEGY_HPP_

#include "Core/Functions.hpp"
#include "Physics/HelicityAmplitude/TwoBodyDecayAmplitude.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

class TwoBodyDecayAngularStrategy: public Strategy {
  std::shared_ptr<TwoBodyDecayAmplitude> tbd_amp_;
  unsigned int storage_index_;

public:
  TwoBodyDecayAngularStrategy(
      std::shared_ptr<TwoBodyDecayAmplitude> tbd_amp,
      unsigned storage_index);
  virtual ~TwoBodyDecayAngularStrategy();

  //! Pure Virtual interface for streaming info about the strategy
  virtual const std::string to_str() const;

  //! Pure Virtual interface for executing a strategy
  virtual bool execute(ParameterList& paras,
      std::shared_ptr<AbsParameter>& out);
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_TWOBODYDECAYANGULARSTRATEGY_HPP_ */
