#ifndef PHYSICS_HELICITYAMPLITUDE_HELICITYAMPLITUDETREEFACTORY_HPP_
#define PHYSICS_HELICITYAMPLITUDE_HELICITYAMPLITUDETREEFACTORY_HPP_

#include "HelicityAmplitudeTree.hpp"

namespace HelicityFormalism {

class HelicityDecayTree;

class HelicityAmplitudeTreeFactory {
public:
  HelicityAmplitudeTreeFactory();
  virtual ~HelicityAmplitudeTreeFactory();

  HelicityAmplitudeTree generateAmplitudeTree(
      const HelicityDecayTree& decay_tree) const;
};

} /* namespace HelicityFormalism */

#endif /* PHYSICS_HELICITYAMPLITUDE_HELICITYAMPLITUDETREEFACTORY_HPP_ */
