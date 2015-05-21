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

#include <algorithm>
#include <stdexcept>

#include "HelicityDecayConfiguration.hpp"

namespace HelicityFormalism {

HelicityDecayConfiguration::HelicityDecayConfiguration() {
}

HelicityDecayConfiguration::~HelicityDecayConfiguration() {
}

void HelicityDecayConfiguration::addFinalStateParticle(
    const ParticleState &final_state_particle) {
  if (std::find(particles_.begin(), particles_.end(), final_state_particle)
      == particles_.end()) {
    particles_.push_back(final_state_particle);
    final_state_particles_.push_back(particles_.size() - 1);
  }
}

void HelicityDecayConfiguration::addIntermediateStateParticle(
    const ParticleState &intermediate_state_particle) {
  if (std::find(particles_.begin(), particles_.end(),
      intermediate_state_particle) == particles_.end()) {
    particles_.push_back(intermediate_state_particle);
    intermediate_state_particles_.push_back(particles_.size() - 1);
  }
}

void HelicityDecayConfiguration::addCurrentDecayTopologyToList() {
  decay_topologies_.push_back(current_decay_topology_);
  current_decay_topology_.clear();
}

unsigned int HelicityDecayConfiguration::convertParticleIDToListIndex(
    unsigned int particle_id) const {
  ParticleStateIDComparator ps_comparator(particle_id);
  std::vector<ParticleState>::const_iterator found_particle = std::find_if(
      particles_.begin(), particles_.end(), ps_comparator);
  if (found_particle != particles_.end()) {
    return std::distance(particles_.begin(), found_particle);
  }
  else {
    throw std::runtime_error(
        "Could not find a particle with the correct ID within the list of particles!");
  }
}

std::vector<unsigned int> HelicityDecayConfiguration::convertParticleIDListToIndexList(
    const std::vector<unsigned int>& particle_id_list) const {
  std::vector<unsigned int> index_list;
  for (unsigned int i = 0; i < particle_id_list.size(); ++i) {
    index_list.push_back(convertParticleIDToListIndex(particle_id_list[i]));
  }
  return index_list;
}

void HelicityDecayConfiguration::addDecayToCurrentDecayTopology(
    unsigned int mother_state_id,
    std::pair<std::vector<unsigned int>, std::vector<unsigned int> > &daughter_states) {
  unsigned int mother_state_index = convertParticleIDToListIndex(mother_state_id);
  current_decay_topology_[mother_state_index].first =
      convertParticleIDListToIndexList(daughter_states.first);
  current_decay_topology_[mother_state_index].second =
      convertParticleIDListToIndexList(daughter_states.second);
}

} /* namespace HelicityFormalism */
