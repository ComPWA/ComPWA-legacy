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

class FullTwoBodyDecayAmplitude : Resonance {
  std::string name;    // for full coherent amplitude construction reasons

  //TODO: we add this particle state info for the coherent sum stuff
  // the whole design is fucked up because of that, change that someday
  std::pair<ParticleStateInfo, std::pair<ParticleStateInfo, ParticleStateInfo> > decay_spin_info_;

  std::shared_ptr<DoubleParameter> strength_;
  std::shared_ptr<DoubleParameter> phase_;
  std::shared_ptr<TwoBodyDecayAmplitude> angular_part_;
  std::shared_ptr<DynamicalFunctions::AbstractDynamicalFunction> dynamical_part_;
  unsigned int positionToEval;

  std::complex<double> Evaluate(dataPoint& point) {
	       sequential_decay_result *= std::polar(
          two_body_decay.strength_->GetValue(),
          two_body_decay.phase_->GetValue());

      std::complex<double> angular_part =
          two_body_decay.angular_part_->evaluate(point,positionToEval);

      std::complex<double> dynamical_part =
          two_body_decay.dynamical_part_->evaluate(point,
              evaluation_index_list[two_body_decay_index]);

      sequential_decay_result *= angular_part * dynamical_part;

  };
};

struct SequentialTwoBodyDecayAmplitude : Amplitude {
  std::vector<FullTwoBodyDecayAmplitude> decay_amplitude;
  double factor;

  SequentialTwoBodyDecayAmplitude() :
      factor(1.0) {
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
