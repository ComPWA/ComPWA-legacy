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

typedef std::pair<std::vector<IDInfo>, IndexList> IDInfoListUniqueIDListPair;

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

  std::vector<ParticleStateInfo> fs_particle_list_;

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

  struct UniqueMappingAccumulatorForTopology {
    const std::set<IndexList>& unique_id_lists_;
    const std::vector<ParticleStateInfo>& fs_particle_list_;
    std::vector<IndexMapping> current_unique_mappings_;

    UniqueMappingAccumulatorForTopology(
        const std::vector<ParticleStateInfo>& fs_particle_list,
        std::set<IndexList>& unique_id_lists) :
        fs_particle_list_(fs_particle_list), unique_id_lists_(unique_id_lists) {
    }

    void addMappingIfUnique(const IndexMapping& mapping) {
      bool add_mapping(true);

      // loop over current mappings
      for (unsigned int i = 0; i < current_unique_mappings_.size(); ++i) {
        // now for each final state content list compare the outcome
        // of the current mapping and the reference
        bool mappings_identical(true);
        for (auto final_state_content_id_list = unique_id_lists_.begin();
            final_state_content_id_list != unique_id_lists_.end();
            ++final_state_content_id_list) {

          if (createFinalStateMappedList(*final_state_content_id_list,
              current_unique_mappings_[i])
              != createFinalStateMappedList(*final_state_content_id_list,
                  mapping)) {
            mappings_identical = false;
            break;
          }
        }
        if (!mappings_identical) {
          break;
        }
        else {
          add_mapping = false;
          break;
        }
      }

      if (add_mapping)
        current_unique_mappings_.push_back(mapping);
    }

    std::set<unsigned int> createFinalStateMappedList(
        const std::vector<unsigned int>& fs_unique_id_list,
        const IndexMapping& mapping) {

      std::set<unsigned int> fs_mapped_list;
      for (unsigned int i = 0; i < fs_unique_id_list.size(); ++i) {
        fs_mapped_list.insert(mapping.at(fs_unique_id_list[i]));
      }

      return fs_mapped_list;
    }
  };

  void createDistinguishableParticleMapping(const Event& event);

  bool isEventParticleMatchUniqueForParticle(const IDInfo& fs_particle,
      const Event& ev) const;

  unsigned int getEventParticleIndex(const IDInfo& particle_state,
      const Event& event) const;

  void createAllParticleMappings(const Event& event);

  void removeDistinguishableParticles(
      std::vector<ParticleStateInfo>& fs_particle_pool,
      std::vector<unsigned int>& event_final_state_particle_index_pool) const;

  void extendParticleMappings(const IndexMapping& current_mapping,
      std::vector<ParticleStateInfo> remaining_final_state_particle_pool,
      const Event& event,
      const std::vector<unsigned int>& remaining_event_final_state_particle_index_pool);

  std::vector<unsigned int> getPossibleEventParticleIndicesForParticleState(
      const IDInfo& particle_state, const Event& event) const;

  void recursiveBuild(unsigned int decay_index,
      const TwoBodyDecayTopology& topology,
      std::map<std::vector<IDInfo>, std::vector<IndexList> >& unique_id_lists,
      const std::vector<ParticleStateInfo>& particle_pool) const;

  std::vector<unsigned int> convertIDInfoToUniqueIDList(
      const std::vector<IDInfo>& fs_content_list,
      std::vector<ParticleStateInfo>& remaining_fs_particles) const;

  std::vector<ParticleStateInfo> createParticleListFromUniqueIDList(
      const std::vector<unsigned int>& unique_id_vector) const;

public:
  FinalStateParticleCombinatorics();
  virtual ~FinalStateParticleCombinatorics();

  void init(const std::vector<ParticleStateInfo>& fs_particle_list, const Event &event);

  std::vector<IndexMapping> getUniqueParticleMappingsSubsetForTopology(
      const TwoBodyDecayTopology& topology) const;

  std::map<std::vector<IDInfo>, std::vector<IndexList> > convertTopologyToUniqueIDLists(
      const TwoBodyDecayTopology& topology) const;
};

} /* namespace HelicityFormalism */

#endif /* PHYSICS_HELICITYAMPLITUDE_FINALSTATEPARTICLECOMBINATORICS_HPP_ */
