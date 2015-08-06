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

#include "TwoBodyDecayAmplitude.hpp"

#include <vector>

namespace HelicityFormalism {

typedef std::vector<
    std::pair<std::shared_ptr<TwoBodyDecayAmplitude>,
        std::shared_ptr<AbstractDynamicalFunction> > > SequentialTwoBodyDecayAmplitude;

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

  std::complex<double> evaluate(
      const std::vector<KinematicVariables>& helicity_angles,
      const IndexList& evaluation_index_list) const;
};

} /* namespace HelicityFormalism */

#endif /* PHYSICS_HELICITYAMPLITUDE_TOPOLOGYAMPLITUDE_HPP_ */
