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

#ifndef PHYSICS_HELICITYAMPLITUDE_FINALSTATEPARTICLECOMBINATORICS_HPP_
#define PHYSICS_HELICITYAMPLITUDE_FINALSTATEPARTICLECOMBINATORICS_HPP_

#include <map>

#include "ParticleStateDefinitions.hpp"

class Event;

namespace HelicityFormalism {

/**
 * This class takes care of the combinatorics that may arise due to
 * indistinguishable final state particles. Its task is to calculate possible
 * mappings of the final state particle states within the amplitude and
 * the event final state particles.
 * The number of mappings depend on the decay topology.
 *
 * The calculation of the mappings is split into two stages (see init()):
 * - At first, all distinguishable final state particles (that can be uniquely
 *   identified) are mapped to the corresponding particle states. This merely
 *   reduces the overall combinatorics.
 * - As a second step, all concrete particle mappings are created. One mapping
 *   for each possible final state particle ordering combination!
 *
 * Since some mappings may be identical to certain decay topologies
 * (i.e: for a f0 decaying into two pi0s, it is irrelevant how the pi0s are
 * assigned since the amplitude is only interested in the center of mass frame
 * information), a subset of unique mappings is obtained for each
 * decay topology via getUniqueParticleMappingsSubsetForTopology().
 */
class FinalStateParticleCombinatorics {
  // one entry for each unique link of a helicity fs particle
  // to an event fs particle
  IndexMapping distinguishable_fs_particle_mapping_;

  std::vector<IndexMapping> fs_particle_mapping_combinations_;

  struct NotInListIndexFilter {
    const std::vector<unsigned int>& list_;

    NotInListIndexFilter(const std::vector<unsigned int>& list) :
        list_(list) {
    }

    bool operator()(unsigned int index) {
      for (unsigned int i = 0; i < list_.size(); ++i) {
        if (index == list_[i])
          return false;
      }
      return true;
    }
  };
  struct InListIndexFilter {
    const std::vector<unsigned int>& list_;

    InListIndexFilter(const std::vector<unsigned int>& list) :
        list_(list) {
    }

    bool operator()(unsigned int index) {
      for (unsigned int i = 0; i < list_.size(); ++i) {
        if (index == list_[i])
          return true;
      }
      return false;
    }
  };
  struct IndexComparison {
    unsigned int index_;

    IndexComparison(unsigned int index) :
        index_(index) {
    }

    bool operator()(unsigned int index) {
      return index_ == index;
    }
  };

  void createDistinguishableParticleMapping(
      const std::vector<ParticleState>& fs_particle_list, const Event& event);

  bool isEventParticleMatchUniqueForParticle(const ParticleState& fs_particle,
      const Event& ev) const;

  unsigned int getEventParticleIndex(const ParticleState& particle_state,
      const Event& event) const;

  void createAllParticleMappings(
      std::vector<ParticleState> final_state_particle_pool, const Event& event);

  void removeDistinguishableParticles(
      std::vector<ParticleState>& fs_particle_pool,
      std::vector<unsigned int>& event_final_state_particle_index_pool) const;

  void extendParticleMappings(
      const IndexMapping& current_mapping,
      std::vector<ParticleState> remaining_final_state_particle_pool,
      const Event& event,
      const std::vector<unsigned int>& remaining_event_final_state_particle_index_pool);

  std::vector<unsigned int> getPossibleEventParticleIndicesForParticleState(
      const ParticleState& particle_state, const Event& event) const;
public:
  FinalStateParticleCombinatorics();
  virtual ~FinalStateParticleCombinatorics();

  void init(const std::vector<ParticleState>& fs_particle_list,
      const Event &event);

  std::vector<IndexMapping> getUniqueParticleMappingsSubsetForTopology(
      const DecayTopology& topology) const;
};

} /* namespace HelicityFormalism */

#endif /* PHYSICS_HELICITYAMPLITUDE_FINALSTATEPARTICLECOMBINATORICS_HPP_ */
