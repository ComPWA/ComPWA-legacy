//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

#ifndef HELICITYDECAYCONFIGURATION_HPP_
#define HELICITYDECAYCONFIGURATION_HPP_

#include "HelicityStateDefinitions.hpp"

#include <map>
#include <vector>

namespace HelicityFormalism {

typedef std::map<unsigned int,
    std::pair<std::vector<unsigned int>, std::vector<unsigned int> > > DecayTopology;

class ParticleStateIDComparator {
  unsigned int particle_id_ ;
public:
  ParticleStateIDComparator(unsigned int particle_id) : particle_id_(particle_id) {
  }
  bool operator()(const ParticleState& ps) const {
    return ps.particle_id_ == particle_id_;
  }
};

class HelicityDecayConfiguration {
  friend class HelicityDecayTreeFactory;

  std::vector<ParticleState> particles_;
  std::vector<unsigned int> final_state_particles_;
  std::vector<unsigned int> intermediate_state_particles_;
  std::vector<DecayTopology> decay_topologies_;
  DecayTopology current_decay_topology_;

public:
  HelicityDecayConfiguration();
  virtual ~HelicityDecayConfiguration();

  void addFinalStateParticle(const ParticleState &final_state_particle);
  void addIntermediateStateParticle(
      const ParticleState &intermediate_state_particle);
  void addCurrentDecayTopologyToList();

  unsigned int convertParticleIDToListIndex(
      unsigned int particle_id) const;

  std::vector<unsigned int> convertParticleIDListToIndexList(
      const std::vector<unsigned int>& particle_id_list) const;

  void addDecayToCurrentDecayTopology(unsigned int mother_state_id,
      std::pair<std::vector<unsigned int>, std::vector<unsigned int> > &daughter_states);
};

} /* namespace HelicityFormalism */

#endif /* HELICITYDECAYCONFIGURATION_HPP_ */
