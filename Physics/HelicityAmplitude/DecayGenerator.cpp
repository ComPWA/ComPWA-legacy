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

// this include has to be after the compwa includes
// because it defines an INTEGER in the global scope that introduces clashes
#include "clips.h"

namespace HelicityFormalism {

DecayGenerator::DecayGenerator() {
  // TODO Auto-generated constructor stub
}

DecayGenerator::~DecayGenerator() {
  // TODO Auto-generated destructor stub
}

void DecayGenerator::addFinalStateParticles(const std::string& name) {
  final_state_particles_.push_back(createIDInfo(name).pid_information_);
  final_state_particles_indices_.push_back(total_particle_pool_.size() - 1);
}
void DecayGenerator::addIntermediateStateParticles(const std::string& name) {
  intermediate_state_particles_.push_back(createIDInfo(name).pid_information_);
  intermediate_state_particles_indices_.push_back(
      total_particle_pool_.size() - 1);
}
void DecayGenerator::setTopNodeState(const std::string& name) {
  mother_state_particle_ = createIDInfo(name).pid_information_;
  mother_state_particle_index_ = total_particle_pool_.size() - 1;
}

ParticleStateInfo DecayGenerator::createIDInfo(const std::string& name) {
  ParticleStateInfo particle_state;
  particle_state.pid_information_.name_ = name;
  particle_state.pid_information_.particle_id_ = lookupPID(name);
  particle_state.unique_id_ = total_particle_pool_.size() + 1;
  total_particle_pool_.push_back(particle_state);
  return particle_state;
}

int DecayGenerator::lookupPID(const std::string& name) const {
  return PhysConst::Instance().findParticle(name).id_;
}

void DecayGenerator::generate() {
  // first construct all possible decay topologies
  /* std::vector<TwoBodyDecayTopology> possible_decay_topologies =
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
   concrete_decay_trees);*/

  void *clips_environment;
  clips_environment = CreateEnvironment();

  //EnvIncrementGCLocks(clips_environment);

  setupClipsEnvironment(clips_environment);

  // run
  int max_number_of_rules_to_fire(-1);
  EnvFacts(clips_environment, "stdout", NULL, -1, -1, -1);
  EnvAgenda(clips_environment, "stdout", NULL);
  EnvRun(clips_environment, max_number_of_rules_to_fire);
  EnvFacts(clips_environment, "stdout", NULL, -1, -1, -1);
  EnvListDefrules(clips_environment, "stdout", NULL);

  // get expertise
  DATA_OBJECT rtn;
  //EnvFunctionCall(theEnv, "getParticle", NULL, &rtn);
  //std::cout << *((char*) rtn.value) << std::endl;

  //cleanup
  //EnvDecrementGCLocks(clips_environment);

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
    trsp.first.insertFinalStateContentList(temp_vec);
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
    std::sort(iter->first.decay_node_fs_content_index_infos_.begin(),
        iter->first.decay_node_fs_content_index_infos_.end());

    if (isTopologyNonExistant(topologies, iter->first)) {
      topologies.push_back(iter->first);
    }
  }

  std::cout << "constructed the " << topologies.size()
      << " decay topologies...\n";
  /*for (auto iter = topologies.begin(); iter != topologies.end(); ++iter) {
   iter->print();
   }*/

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
    for (unsigned int index_remaining_state = 0;
        index_remaining_state < topology_remaining_state_pair->second.size();
        ++index_remaining_state) {
      for (unsigned int index_second_remaining_state = index_remaining_state
          + 1;
          index_second_remaining_state
              < topology_remaining_state_pair->second.size();
          ++index_second_remaining_state) {
        //create a new TopologyRemainingStatePair
        TopologyRemainingStatePair trsp(*topology_remaining_state_pair);

        TwoBodyDecayIndices indices;

        indices.decay_products_.first = trsp.first.findFinalStateList(
            trsp.second[index_remaining_state]);

        indices.decay_products_.second = trsp.first.findFinalStateList(
            trsp.second[index_second_remaining_state]);

        trsp.second[index_remaining_state].insert(
            trsp.second[index_remaining_state].end(),
            trsp.second[index_second_remaining_state].begin(),
            trsp.second[index_second_remaining_state].end());

        indices.decay_state_index_ = trsp.first.insertFinalStateContentList(
            trsp.second[index_remaining_state]);

        trsp.first.setMotherDecayIndexForDecayIndex(indices.decay_state_index_,
            indices.decay_products_.first);
        trsp.first.setMotherDecayIndexForDecayIndex(indices.decay_state_index_,
            indices.decay_products_.second);

        // we set the mother state as the decay state because mother cannot
        // exist yet and will be linked later correctly
        indices.mother_index_ = indices.decay_state_index_;
        //trsp.second.erase(selected_state_1);
        trsp.second.erase(trsp.second.begin() + index_second_remaining_state);

        trsp.first.decay_node_fs_content_index_infos_.push_back(indices);

        new_topology_remaining_state_pairs.push_back(trsp);
      }
    }
  }

  return new_topology_remaining_state_pairs;
}
/*
std::vector<ParticleIndexDecayTree> DecayGenerator::createDecayTreeRealizationTemplates(
    const std::vector<TwoBodyDecayTopology>& possible_decay_topologies) const {
  std::vector<ParticleIndexDecayTree> decay_tree_realization_templates;

  return decay_tree_realization_templates;
}

std::vector<TParticlePDG> DecayGenerator::createDecayNodeParticleCandiateList(
    const std::vector<IDInfo>& decay_products) const {
}

std::vector<TParticlePDG> DecayGenerator::convertIDInfoToTParticlePDG(
    const std::vector<IDInfo>& particle_id_info_list) const {

  for (auto particle_id_info = particle_id_info_list.begin();
   particle_id_info != particle_id_info_list.end(); ++particle_id_info) {
   particle_list.push_back(TParticlePDG()particle_id_info->particle_id_));
   }
  return particle_list;
}*/

void DecayGenerator::setupClipsEnvironment(void *clips_environment) const {
  std::stringstream path;
  path << getenv("COMPWA_DIR")
      << "/Physics/HelicityAmplitude/ConcreteDecayRules/helicity_model.clp";
  EnvLoad(clips_environment, path.str().c_str());
  //EnvReset(clips_environment);

  for (unsigned int i = 0; i < total_particle_pool_.size(); ++i) {
    addParticleToClipsEnvironment(clips_environment, total_particle_pool_[i]);
  }

  // create top node global unique id
  DATA_OBJECT top_node_id_global_value;
  top_node_id_global_value.type = INTEGER;
  top_node_id_global_value.value = EnvAddLong(clips_environment,
      total_particle_pool_[mother_state_particle_index_].unique_id_);
  EnvSetDefglobalValue(clips_environment, "initial_state_unique_id",
      &top_node_id_global_value);

  // create two body decay tree seed
  void* twobodydecaytree_template;
  twobodydecaytree_template = EnvFindDeftemplate(clips_environment,
      "TwoBodyDecayTree");
  void* twobodydecaytree_fact = EnvCreateFact(clips_environment,
      twobodydecaytree_template);
  if (twobodydecaytree_fact != NULL) {
    DATA_OBJECT available_final_state_particles;
    void* multifield_ptr = EnvCreateMultifield(clips_environment,
        final_state_particles_.size());

    for (unsigned int i = 0; i < final_state_particles_indices_.size(); ++i) {
      SetMFType(multifield_ptr, i + 1, INTEGER);
      SetMFValue(multifield_ptr, i + 1,
          EnvAddLong(clips_environment,
              total_particle_pool_[final_state_particles_indices_[i]].unique_id_));
    }

    SetpType(&available_final_state_particles, MULTIFIELD);
    SetpValue(&available_final_state_particles, multifield_ptr);

    SetpDOBegin(&available_final_state_particles, 1);
    SetpDOEnd(&available_final_state_particles,
        final_state_particles_indices_.size());

    EnvPutFactSlot(clips_environment, twobodydecaytree_fact,
        "available_particles", &available_final_state_particles);

   /* DATA_OBJECT two_body_decays;
    multifield_ptr = EnvCreateMultifield(clips_environment, 0);
    SetpType(&two_body_decays, MULTIFIELD);
    SetpValue(&two_body_decays, multifield_ptr);

    SetpDOBegin(&two_body_decays, 0);
    SetpDOEnd(&two_body_decays, 0);
    EnvPutFactSlot(clips_environment, twobodydecaytree_fact, "two_body_decays",
        &two_body_decays);*/

    EnvAssignFactSlotDefaults(clips_environment, twobodydecaytree_fact);

    EnvAssert(clips_environment, twobodydecaytree_fact);
  }

  // create result fact
  void* result_template;
  result_template = EnvFindDeftemplate(clips_environment,
      "ValidTwoBodyDecayTrees");
  void* result_fact = EnvCreateFact(clips_environment,
      result_template);
  EnvAssignFactSlotDefaults(clips_environment, result_fact);
  EnvAssert(clips_environment, result_fact);
}

/* a particle in clips consists of the following slots
 (slot unique_id (type INTEGER))
 (slot name (type STRING))
 (slot pid (type INTEGER))
 (slot mass (type FLOAT))
 (slot charge (type INTEGER))
 (slot isospin (type INTEGER))
 (slot spin (type INTEGER))
 (slot parity (type INTEGER))
 (slot cparity (type INTEGER))
 */
void DecayGenerator::addParticleToClipsEnvironment(void *clips_environment,
    const ParticleStateInfo& particle_info) const {
  ParticleProperties particle_properties = PhysConst::Instance().findParticle(
      particle_info.pid_information_.particle_id_);

  void* particle_template;
  particle_template = EnvFindDeftemplate(clips_environment, "Particle");
  // set the facts
  void* particle_fact = EnvCreateFact(clips_environment, particle_template);
  if (particle_fact != NULL) {
    DATA_OBJECT data_value;
    data_value.type = INTEGER;
    data_value.value = EnvAddLong(clips_environment, particle_info.unique_id_);
    EnvPutFactSlot(clips_environment, particle_fact, "unique_id", &data_value);
    data_value.type = SYMBOL;
    data_value.value = EnvAddSymbol(clips_environment,
        particle_properties.name_.c_str());
    EnvPutFactSlot(clips_environment, particle_fact, "name", &data_value);
    data_value.type = INTEGER;
    data_value.value = EnvAddLong(clips_environment, particle_properties.id_);
    EnvPutFactSlot(clips_environment, particle_fact, "pid", &data_value);
    data_value.type = FLOAT;
    data_value.value = EnvAddDouble(clips_environment,
        particle_properties.mass_);
    EnvPutFactSlot(clips_environment, particle_fact, "mass", &data_value);
    data_value.type = INTEGER;
    data_value.value = EnvAddLong(clips_environment,
        particle_properties.charge_);
    EnvPutFactSlot(clips_environment, particle_fact, "charge", &data_value);
    data_value.value = EnvAddLong(clips_environment,
        particle_properties.isospin_);
    EnvPutFactSlot(clips_environment, particle_fact, "isospin", &data_value);
    data_value.value = EnvAddLong(clips_environment,
        particle_properties.isospin_z_);
    EnvPutFactSlot(clips_environment, particle_fact, "isospin_z", &data_value);
    data_value.value = EnvAddLong(clips_environment, particle_properties.spin_);
    EnvPutFactSlot(clips_environment, particle_fact, "spin", &data_value);
    data_value.value = EnvAddLong(clips_environment,
        particle_properties.parity_);
    EnvPutFactSlot(clips_environment, particle_fact, "parity", &data_value);
    data_value.value = EnvAddLong(clips_environment,
        particle_properties.cparity_);
    EnvPutFactSlot(clips_environment, particle_fact, "cparity", &data_value);
  }

  EnvAssert(clips_environment, particle_fact);
}

std::vector<ParticleIndexDecayTree> DecayGenerator::makeConcreteDecayTrees(
    const std::vector<ParticleIndexDecayTree>& decay_tree_realization_templates) {
  std::vector<ParticleIndexDecayTree> concrete_decay_trees;

  return concrete_decay_trees;
}

boost::property_tree::ptree DecayGenerator::portPropertyTree(
    const std::vector<ParticleIndexDecayTree>& concrete_decay_trees) const {
  return boost::property_tree::ptree();
}

} /* namespace HelicityFormalism */
