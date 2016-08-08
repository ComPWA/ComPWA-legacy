//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//----------------------------------------------------------------------------------

#ifndef PHYSICS_HELICITYAMPLITUDE_TOPOLOGYAMPLITUDE_HPP_
#define PHYSICS_HELICITYAMPLITUDE_TOPOLOGYAMPLITUDE_HPP_

#include <vector>

#include "TwoBodyDecayAmplitude.hpp"

#include "Physics/DynamicalDecayFunctions/AbstractDynamicalFunction.hpp"

namespace ComPWA {

class dataPoint;

namespace Physics {
namespace HelicityFormalism {

struct FullTwoBodyDecayAmplitude {
  std::shared_ptr<DoubleParameter> strength_;
  std::shared_ptr<DoubleParameter> phase_;
  std::shared_ptr<TwoBodyDecayAmplitude> angular_part_;
  std::shared_ptr<DynamicalFunctions::AbstractDynamicalFunction> dynamical_part_;
};

struct SequentialTwoBodyDecayAmplitude {
  std::vector<FullTwoBodyDecayAmplitude> decay_amplitude;
  double factor;

  SequentialTwoBodyDecayAmplitude() : factor(1.0) {
  }
};

/**
 * This class defines the full amplitude constructed of all decay amplitudes
 * of the same decay topology. They all share the same kinematics, the helicity
 * angles in each decay node.
 */
class TopologyAmplitude {
  friend class TopologyAmplitudeFactory;

  std::vector<SequentialTwoBodyDecayAmplitude> sequential_decay_amplitude_list_;

public:
  TopologyAmplitude();
  virtual ~TopologyAmplitude();

  std::complex<double> evaluate(const dataPoint& point,
      const IndexList& evaluation_index_list) const;

  const std::vector<SequentialTwoBodyDecayAmplitude>& getSequentialDecayList() const;
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_TOPOLOGYAMPLITUDE_HPP_ */
