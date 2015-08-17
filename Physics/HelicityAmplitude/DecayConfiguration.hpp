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

#ifndef PHYSICS_HELICITYAMPLITUDE_DECAYCONFIGURATION_HPP_
#define PHYSICS_HELICITYAMPLITUDE_DECAYCONFIGURATION_HPP_

#include "ParticleStateDefinitions.hpp"

#include <map>
#include <vector>

namespace HelicityFormalism {

typedef std::map<unsigned int, std::vector<std::vector<unsigned int> > > ParticleIndexDecayTree;

class ParticleStateIDComparator {
  unsigned int particle_id_;
public:
  ParticleStateIDComparator(unsigned int particle_id) :
      particle_id_(particle_id) {
  }
  bool operator()(const ParticleStateInfo& ps) const {
    return ps.id_ == particle_id_;
  }
};

class DecayConfiguration {
  friend class DecayTreeFactory;

  std::vector<ParticleStateInfo> particles_;
  std::vector<unsigned int> final_state_particles_;
  std::vector<unsigned int> intermediate_state_particles_;

  std::vector<ParticleIndexDecayTree> concrete_decay_trees_;

  ParticleIndexDecayTree current_concrete_decay_tree_;

public:
  DecayConfiguration();
  virtual ~DecayConfiguration();

  void addFinalStateParticle(const ParticleStateInfo &final_state_particle);
  void addIntermediateStateParticle(
      const ParticleStateInfo &intermediate_state_particle);
  void addCurrentDecayTreeToList();

  unsigned int convertParticleIDToListIndex(unsigned int particle_id) const;

  std::vector<unsigned int> convertParticleIDListToIndexList(
      const std::vector<unsigned int>& particle_id_list) const;

  void addDecayToCurrentDecayTree(unsigned int mother_state_id,
      std::vector<std::vector<unsigned int> > &daughter_states);
};

} /* namespace HelicityFormalism */

#endif /* PHYSICS_HELICITYAMPLITUDE_DECAYCONFIGURATION_HPP_ */
