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

#include "FinalStateParticleCombinatorics.hpp"

#include "Core/Event.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

FinalStateParticleCombinatorics::FinalStateParticleCombinatorics() {
}

FinalStateParticleCombinatorics::~FinalStateParticleCombinatorics() {
}

void FinalStateParticleCombinatorics::init(
    const std::vector<ParticleStateInfo>& fs_particle_list,
    const Event& event) {
  fs_particle_list_ = fs_particle_list;

  createDistinguishableParticleMapping(event);

  createAllParticleMappings(event);
}

void FinalStateParticleCombinatorics::createDistinguishableParticleMapping(
    const Event& event) {
  for (auto fs_particle_iter = fs_particle_list_.begin();
      fs_particle_iter != fs_particle_list_.end(); ++fs_particle_iter) {
    if (isEventParticleMatchUniqueForParticle(
        fs_particle_iter->pid_information_, event)) {
      distinguishable_fs_particle_mapping_[fs_particle_iter->unique_id_] =
          getEventParticleIndex(fs_particle_iter->pid_information_, event);
    }
  }
}

bool FinalStateParticleCombinatorics::isEventParticleMatchUniqueForParticle(
    const IDInfo& fs_particle, const Event& event) const {
  unsigned int match_counter(0);
  for (unsigned int i = 0; i < event.getNParticles(); ++i) {
    if (event.getParticle(i).pid == fs_particle.particle_id_)
      ++match_counter;
  }

  if (match_counter == 1)
    return true;
  else
    return false;
}

unsigned int FinalStateParticleCombinatorics::getEventParticleIndex(
    const IDInfo& particle_state, const Event& event) const {
  for (unsigned int i = 0; i < event.getNParticles(); ++i) {
    if (event.getParticle(i).pid == particle_state.particle_id_)
      return i;
  }

  std::stringstream ss;
  ss
      << "FinalStateParticleCombinatorics: unable to find a matching particle in the event structure for "
      << particle_state.name_ << " with PID " << particle_state.particle_id_;
  throw std::runtime_error(ss.str());
}

void FinalStateParticleCombinatorics::createAllParticleMappings(
    const Event& event) {

  std::vector<ParticleStateInfo> final_state_particle_pool(fs_particle_list_);

  // create event particle index pool
  std::vector<unsigned int> event_final_state_particle_index_pool;
  for (unsigned int i = 0; i < event.getNParticles(); ++i) {
    event_final_state_particle_index_pool.push_back(i);
  }
  // remove the distinguishable particle mapping particles
  removeDistinguishableParticles(final_state_particle_pool,
      event_final_state_particle_index_pool);

  // initialize the starting mapping with the distinguishable particle mapping
  IndexMapping start_mapping(distinguishable_fs_particle_mapping_);

  // then start the recursive filling algorithm with these filtered particle pools
  extendParticleMappings(start_mapping, final_state_particle_pool, event,
      event_final_state_particle_index_pool);
}

void FinalStateParticleCombinatorics::removeDistinguishableParticles(
    std::vector<ParticleStateInfo>& fs_particle_pool,
    std::vector<unsigned int>& event_final_state_particle_index_pool) const {

  IndexMapping::const_iterator mapping_link_iter;
  for (mapping_link_iter = distinguishable_fs_particle_mapping_.begin();
      mapping_link_iter != distinguishable_fs_particle_mapping_.end();
      ++mapping_link_iter) {
    fs_particle_pool.erase(
        std::remove_if(fs_particle_pool.begin(), fs_particle_pool.end(),
            ParticleIDComparison(mapping_link_iter->first)),
        fs_particle_pool.end());
    event_final_state_particle_index_pool.erase(
        std::remove_if(event_final_state_particle_index_pool.begin(),
            event_final_state_particle_index_pool.end(),
            IndexComparison(mapping_link_iter->second)),
        event_final_state_particle_index_pool.end());
  }
}

void FinalStateParticleCombinatorics::extendParticleMappings(
    const IndexMapping& current_mapping,
    std::vector<ParticleStateInfo> remaining_final_state_particle_pool,
    const Event& event,
    const std::vector<unsigned int>& remaining_event_final_state_particle_index_pool) {

  if (remaining_final_state_particle_pool.size() > 0) {
    const ParticleStateInfo &ps = remaining_final_state_particle_pool.back();
    remaining_final_state_particle_pool.pop_back();

    // first get all candidates for the particle
    std::vector<unsigned int> event_fs_candidates =
        getPossibleEventParticleIndicesForParticleState(ps.pid_information_,
            event);

    // filter only for possible remaining candidates
    event_fs_candidates.erase(
        std::remove_if(event_fs_candidates.begin(), event_fs_candidates.end(),
            NotInListIndexFilter(
                remaining_event_final_state_particle_index_pool)),
        event_fs_candidates.end());

    //extend the current mappings
    std::vector<std::pair<IndexMapping, std::vector<unsigned int> > > extended_mapping_stripped_pool_pairs;

    for (unsigned int i = 0; i < event_fs_candidates.size(); ++i) {
      IndexMapping extended_mapping(current_mapping);

      extended_mapping[ps.unique_id_] = event_fs_candidates[i];

      std::vector<unsigned int> temp_vector;
      temp_vector.push_back(event_fs_candidates[i]);
      std::vector<unsigned int> stripped_pool(
          remaining_event_final_state_particle_index_pool);
      stripped_pool.erase(
          std::remove_if(stripped_pool.begin(), stripped_pool.end(),
              InListIndexFilter(temp_vector)), stripped_pool.end());

      extended_mapping_stripped_pool_pairs.push_back(
          std::make_pair(extended_mapping, stripped_pool));
    }

    // and recurse
    for (unsigned int i = 0; i < extended_mapping_stripped_pool_pairs.size();
        ++i) {
      extendParticleMappings(extended_mapping_stripped_pool_pairs[i].first,
          remaining_final_state_particle_pool, event,
          extended_mapping_stripped_pool_pairs[i].second);
    }
  }
  else {
    fs_particle_mapping_combinations_.insert(
        fs_particle_mapping_combinations_.begin(), current_mapping);
  }
}

std::vector<unsigned int> FinalStateParticleCombinatorics::getPossibleEventParticleIndicesForParticleState(
    const IDInfo& particle_id_info, const Event& event) const {
  std::vector<unsigned int> index_list;
  for (unsigned int i = 0; i < event.getNParticles(); ++i) {
    if (event.getParticle(i).pid == particle_id_info.particle_id_)
      index_list.push_back(i);
  }
  return index_list;
}

std::vector<IndexMapping> FinalStateParticleCombinatorics::getUniqueParticleMappingsSubsetForTopology(
    const TwoBodyDecayTopology& topology) const {

  //convert the fs particle lists from the topology to
  //a list of particle unique ids
  std::set<IndexList> unique_id_lists;
  std::map<unsigned int, IndexList> unique_id_map =
      topology.final_state_content_unique_id_mapping_;
  for (auto iter = unique_id_map.begin(); iter != unique_id_map.end(); ++iter) {
    if (iter->second.size() > 1)
      unique_id_lists.insert(iter->second);
  }

  UniqueMappingAccumulatorForTopology mapping_accumulator(fs_particle_list_,
      unique_id_lists);

  std::vector<IndexMapping>::const_iterator index_mapping_iter;

  for (index_mapping_iter = fs_particle_mapping_combinations_.begin();
      index_mapping_iter != fs_particle_mapping_combinations_.end();
      ++index_mapping_iter) {

    mapping_accumulator.addMappingIfUnique(*index_mapping_iter);
  }
  return mapping_accumulator.current_unique_mappings_;
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
