//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//    Stefan Pflueger - initial implementation
//-------------------------------------------------------------------------------

#include "Core/PhysConst.hpp"

#include "Physics/HelicityAmplitude/DecayGenerator.hpp"

namespace HelicityFormalism {

DecayGenerator::DecayGenerator() {
  // TODO Auto-generated constructor stub

}

DecayGenerator::~DecayGenerator() {
  // TODO Auto-generated destructor stub
}

void DecayGenerator::addFinalStateParticles(const std::string& name) {
  final_state_particles_.push_back(createIDInfo(name));
}
void DecayGenerator::addIntermediateStateParticles(const std::string& name) {
  intermediate_state_particles_.push_back(createIDInfo(name));
}
void DecayGenerator::setTopNodeState(const std::string& name) {
  mother_state_particle_ = createIDInfo(name);
}

IDInfo DecayGenerator::createIDInfo(const std::string& name) {
  IDInfo pinfo;
  pinfo.name_ = name;
  pinfo.particle_id_ = lookupPID(name);
  pinfo.id_ = total_particle_pool_.size() + 1;
  total_particle_pool_.push_back(pinfo);
  return pinfo;
}

int DecayGenerator::lookupPID(const std::string& name) const {
  PhysConst *physics_constants = PhysConst::instance();
  return physics_constants->getId(name);
}

void DecayGenerator::generate() {
  // first construct all possible decay topologies
  std::vector<TwoBodyDecayTopology> possible_decay_topologies =
      constructPossibleDecayTopologies();

  // then check the intermediate state pool for possible realization templates
  // of the different decay topologies
  std::vector<ParticleIndexDecayTree> decay_tree_realization_templates =
      createDecayTreeRealizationTemplates(possible_decay_topologies);

  // for each realization template, just generate concrete realizations with
  // all allowed M L combinations based on the user input
  std::vector<ParticleIndexDecayTree> concrete_decay_trees =
      makeConcreteDecayTrees(decay_tree_realization_templates);

  // now just port all the concrete realizations to an boost property tree
  // and write them to a xml file
  boost::property_tree::ptree decay_property_tree = portPropertyTree(
      concrete_decay_trees);

  // thats it!
}

std::vector<TwoBodyDecayTopology> DecayGenerator::constructPossibleDecayTopologies() const {
  // initialize
  std::vector<TopologyRemainingStatePair> current_topology_remaining_state_pairs;
  TopologyRemainingStatePair trsp;
  trsp.first.final_state_id_list_ = final_state_particles_;
  trsp.first.top_node_id_info_ = mother_state_particle_;
  for (auto fs_particle = final_state_particles_.begin();
      fs_particle != final_state_particles_.end(); ++fs_particle) {
    std::vector<IDInfo> temp_vec;
    temp_vec.push_back(*fs_particle);
    trsp.second.push_back(temp_vec);
  }
  current_topology_remaining_state_pairs.push_back(trsp);

  while (!isTopologyBuildingFinished(current_topology_remaining_state_pairs)) {
    current_topology_remaining_state_pairs = constructNextLevel(
        current_topology_remaining_state_pairs);
  }

  std::vector<TwoBodyDecayTopology> topologies;
  for (auto iter = current_topology_remaining_state_pairs.begin();
      iter != current_topology_remaining_state_pairs.end(); ++iter) {
    std::sort(iter->first.final_state_content_id_lists_.begin(),
        iter->first.final_state_content_id_lists_.end());
    std::sort(iter->first.decay_node_infos_.begin(),
        iter->first.decay_node_infos_.end());
    if (isTopologyNonExistant(topologies, iter->first)) {
      topologies.push_back(iter->first);
    }
  }

  std::cout << "constructed the following decay topologies...\n";
  for (auto iter = topologies.begin(); iter != topologies.end(); ++iter) {
    iter->print();
  }

  return topologies;
}

bool DecayGenerator::isTopologyBuildingFinished(
    const std::vector<TopologyRemainingStatePair>& current_topology_remaining_state_pairs) const {
  for (auto iter = current_topology_remaining_state_pairs.begin();
      iter != current_topology_remaining_state_pairs.end(); ++iter) {
    if (iter->second.size() > 1)
      return false;
  }

  return true;
}

bool DecayGenerator::isTopologyNonExistant(
    const std::vector<TwoBodyDecayTopology> topology_pool,
    const TwoBodyDecayTopology& probe) const {
  if (std::find(topology_pool.begin(), topology_pool.end(), probe)
      != topology_pool.end()) {
    return false;
  }
  return true;
}

std::vector<TopologyRemainingStatePair> DecayGenerator::constructNextLevel(
    const std::vector<TopologyRemainingStatePair>& current_topology_remaining_state_pairs) const {

  std::vector<TopologyRemainingStatePair> new_topology_remaining_state_pairs;

  for (auto topology_remaining_state_pair =
      current_topology_remaining_state_pairs.begin();
      topology_remaining_state_pair
          != current_topology_remaining_state_pairs.end();
      ++topology_remaining_state_pair) {

    // for each topology make all possible remaining state groupings and
    for (auto remaining_state = topology_remaining_state_pair->second.begin();
        remaining_state != topology_remaining_state_pair->second.end();
        ++remaining_state) {
      auto second_remaining_state = remaining_state;
      for (std::advance(second_remaining_state, 1);
          second_remaining_state != topology_remaining_state_pair->second.end();
          ++second_remaining_state) {
        //create a new TopologyRemainingStatePair
        TopologyRemainingStatePair trsp(*topology_remaining_state_pair);

        std::vector<std::vector<IDInfo> >::iterator selected_state_1 =
            std::find(trsp.second.begin(), trsp.second.end(), *remaining_state);
        std::vector<std::vector<IDInfo> >::iterator selected_state_2 =
            std::find(trsp.second.begin(), trsp.second.end(),
                *second_remaining_state);

        if (selected_state_1 != trsp.second.end()
            && selected_state_2 != trsp.second.end()) {
          std::vector<IDInfo> new_fs_particle_list(*selected_state_1);
          new_fs_particle_list.insert(new_fs_particle_list.end(),
              selected_state_2->begin(), selected_state_2->end());

          trsp.first.insertFinalStateContentList(new_fs_particle_list);

          trsp.second.erase(selected_state_1);
          trsp.second.erase(selected_state_2);
          trsp.second.push_back(new_fs_particle_list);

          new_topology_remaining_state_pairs.push_back(trsp);
        }
      }
    }
  }

  return new_topology_remaining_state_pairs;
}

std::vector<ParticleIndexDecayTree> DecayGenerator::createDecayTreeRealizationTemplates(
    const std::vector<TwoBodyDecayTopology>& possible_decay_topologies) const {
  std::vector<ParticleIndexDecayTree> decay_tree_realization_templates;

  return decay_tree_realization_templates;
}

std::vector<ParticleIndexDecayTree> DecayGenerator::makeConcreteDecayTrees(
    const std::vector<ParticleIndexDecayTree>& decay_tree_realization_templates) {
  std::vector<ParticleIndexDecayTree> concrete_decay_trees;

  return concrete_decay_trees;
}

boost::property_tree::ptree DecayGenerator::portPropertyTree(
    const std::vector<ParticleIndexDecayTree>& concrete_decay_trees) const {

}

} /* namespace HelicityFormalism */
