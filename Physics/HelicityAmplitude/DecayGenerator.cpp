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
  clips_environment = CreateEnvironment();
}

DecayGenerator::~DecayGenerator() {
  // TODO Auto-generated destructor stub
}

IFParticleInfo DecayGenerator::createIFParticleInfo(
    const std::string& name) const {
  IFParticleInfo particle_state;
  particle_state.particle_info_.name_ = name;
  particle_state.particle_info_.particle_id_ = lookupPID(name);
  return particle_state;
}

void DecayGenerator::addFinalStateParticles(
    const IFParticleInfo& particle_info) {
  final_state_particles_.push_back(particle_info);
  IFParticleInfo& pi = final_state_particles_[final_state_particles_.size() - 1];
  pi.unique_id_ = final_state_particles_.size();    // 0 is reserved for the top node
  populateSpinZStates(pi);
  //IndexList indices = createParticleStates(pi);
  //final_state_particles_indices_.insert(final_state_particles_indices_.end(),
  //    indices.begin(), indices.end());
}

/*void DecayGenerator::addIntermediateStateParticles(const std::string& name) {
 intermediate_state_particles_.push_back(createIDInfo(name).pid_information_);
 intermediate_state_particles_indices_.push_back(
 total_particle_pool_.size() - 1);
 }*/

void DecayGenerator::setTopNodeState(const IFParticleInfo& particle_info) {
  mother_state_particle_ = particle_info;
  mother_state_particle_.unique_id_ = 0;
  populateSpinZStates(mother_state_particle_);
  //IndexList indices = createParticleStates(mother_state_particle_);
  //mother_state_particle_index_.insert(mother_state_particle_index_.end(),
  //    indices.begin(), indices.end());

}

void DecayGenerator::populateSpinZStates(IFParticleInfo& particle_info) const {
  // if the user did not specify any spin z components yet
  // then just populate homogeneously
  if (particle_info.spin_z_components_.size() == 0) {
    ComPWA::Spin s = ComPWA::PhysConst::Instance().findParticle(
        particle_info.particle_info_.name_).spin_;

    for (int i = -1 * s.J_numerator_; i <= (int) s.J_numerator_;
        i = i + s.J_denominator_) {
      s.J_z_numerator_ = i;
      particle_info.spin_z_components_.push_back(s);
    }
  }
}

int DecayGenerator::lookupPID(const std::string& name) const {
  return ComPWA::PhysConst::Instance().findParticle(name).id_;
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
   */


  createAllSpinWaves();
  createAllTwoBodyDecays();
  validateTwoBodyDecays();
  createAllValidTwoBodyDecayTrees();

  //cleanup clips environment
  //  EnvDecrementGCLocks(clips_environment);

  std::vector<ParticleIndexDecayTree> concrete_decay_trees;
  // now just port all the concrete realizations to an boost property tree
  // and write them to a xml file
  boost::property_tree::ptree decay_property_tree = portPropertyTree(
      concrete_decay_trees);
}

void DecayGenerator::createAllSpinWaves() {

}

void DecayGenerator::createAllTwoBodyDecays() {

}

void DecayGenerator::validateTwoBodyDecays() {
  //  EnvIncrementGCLocks(clips_environment);
  setupClipsEnvironment();

  // set other user conditions (allowed quantum numbers)
  applyUserConditions();

  // run
  //EnvFacts(clips_environment, "stdout", NULL, -1, -1, -1);
  int max_number_of_rules_to_fire(-1);
  EnvRun(clips_environment, max_number_of_rules_to_fire);
  EnvFacts(clips_environment, "stdout", NULL, -1, -1, -1);

  // get expertise
  getExpertise();
}

/*
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
/*
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

void DecayGenerator::setupClipsEnvironment() const {
  std::stringstream path;
  path << getenv("COMPWA_DIR")
      << "/Physics/HelicityAmplitude/helicity_model.clp";
  EnvLoad(clips_environment, path.str().c_str());
  EnvReset(clips_environment);

  addIFParticleToClipsEnvironment(mother_state_particle_);

  for (unsigned int i = 0; i < final_state_particles_.size(); ++i) {
    addIFParticleToClipsEnvironment(final_state_particles_[i]);
  }

// create top node global unique id
//  DATA_OBJECT top_node_id_global_value;
//  top_node_id_global_value.type = INTEGER;
//  top_node_id_global_value.value = EnvAddLong(clips_environment,
//      total_particle_pool_[mother_state_particle_index_].unique_id_);
//  EnvSetDefglobalValue(clips_environment, "initial_state_unique_id",
//      &top_node_id_global_value);

// create two body decay tree seed
  /*  void* twobodydecaytree_template;
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

  /* EnvAssignFactSlotDefaults(clips_environment, twobodydecaytree_fact);

   EnvAssert(clips_environment, twobodydecaytree_fact);
   }*/
}

void DecayGenerator::addIFParticleToClipsEnvironment(
    const IFParticleInfo& particle_info) const {
  ComPWA::ParticleProperties particle_properties =
      ComPWA::PhysConst::Instance().findParticle(
          particle_info.particle_info_.name_);

  void* particle_template;
  particle_template = EnvFindDeftemplate(clips_environment,
      "SpinWaveMultiplet");
// set the facts
  void* particle_fact = EnvCreateFact(clips_environment, particle_template);
  if (particle_fact != NULL) {
    DATA_OBJECT data_value;
    data_value.type = INTEGER;
    data_value.value = EnvAddLong(clips_environment, particle_info.unique_id_);
    EnvPutFactSlot(clips_environment, particle_fact, "unique_id", &data_value);
    /*data_value.type = SYMBOL;
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
     data_value.type = INTEGER;*/
    data_value.value = EnvAddLong(clips_environment,
        particle_properties.charge_);
    EnvPutFactSlot(clips_environment, particle_fact, "charge", &data_value);
    data_value.value = EnvAddLong(clips_environment,
        particle_properties.isospin_.J_numerator_);
    EnvPutFactSlot(clips_environment, particle_fact, "isospin_num",
        &data_value);
    data_value.value = EnvAddLong(clips_environment,
        particle_properties.isospin_.J_denominator_);
    EnvPutFactSlot(clips_environment, particle_fact, "isospin_denom",
        &data_value);
    data_value.value = EnvAddLong(clips_environment,
        particle_properties.isospin_.J_z_numerator_);
    EnvPutFactSlot(clips_environment, particle_fact, "isospin_z_num",
        &data_value);
    data_value.value = EnvAddLong(clips_environment,
        particle_properties.spin_.J_numerator_);
    EnvPutFactSlot(clips_environment, particle_fact, "spin_num", &data_value);
    data_value.value = EnvAddLong(clips_environment,
        particle_properties.spin_.J_denominator_);
    EnvPutFactSlot(clips_environment, particle_fact, "spin_denom", &data_value);

    DATA_OBJECT spin_z_components;
    void* multifield_ptr = EnvCreateMultifield(clips_environment,
        particle_info.spin_z_components_.size());

    for (unsigned int i = 0; i < particle_info.spin_z_components_.size(); ++i) {
      SetMFType(multifield_ptr, i + 1, INTEGER);
      SetMFValue(multifield_ptr, i + 1,
          EnvAddLong(clips_environment,
              particle_info.spin_z_components_[i].J_z_numerator_));
    }

    SetType(spin_z_components, MULTIFIELD);
    SetValue(spin_z_components, multifield_ptr);

    SetDOBegin(spin_z_components, 1);
    SetDOEnd(spin_z_components, particle_info.spin_z_components_.size());

    EnvPutFactSlot(clips_environment, particle_fact, "spin_z_num",
        &spin_z_components);

    data_value.value = EnvAddLong(clips_environment,
        particle_properties.parity_);
    EnvPutFactSlot(clips_environment, particle_fact, "parity", &data_value);
    data_value.value = EnvAddLong(clips_environment,
        particle_properties.cparity_);
    EnvPutFactSlot(clips_environment, particle_fact, "cparity", &data_value);

    EnvAssignFactSlotDefaults(clips_environment, particle_fact);
    EnvAssert(clips_environment, particle_fact);
  }
}

void DecayGenerator::applyUserConditions() const {
  /*(AllowedQN
        (spin_nums 0 1 2) (spin_denom 1 ) (isospin_nums 0 1) (isospin_denom 1)
        (charge -1 0 1) (parity -1 1) (cparity -1 1)
    )

  void* allowedqn_template;
  allowedqn_template = EnvFindDeftemplate(clips_environment,
      "AllowedQN");
  void* twobodydecaytree_fact = EnvCreateFact(clips_environment,
      allowedqn_template);
  if (twobodydecaytree_fact != NULL) {
    DATA_OBJECT available_final_state_particles;
    void* multifield_ptr = EnvCreateMultifield(clips_environment,
        allowed_charges_.size());

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

    DATA_OBJECT two_body_decays;
    multifield_ptr = EnvCreateMultifield(clips_environment, 0);
    SetpType(&two_body_decays, MULTIFIELD);
    SetpValue(&two_body_decays, multifield_ptr);

    SetpDOBegin(&two_body_decays, 0);
    SetpDOEnd(&two_body_decays, 0);
    EnvPutFactSlot(clips_environment, twobodydecaytree_fact, "two_body_decays",
        &two_body_decays);

    EnvAssignFactSlotDefaults(clips_environment, twobodydecaytree_fact);
  }*/
}

void DecayGenerator::getExpertise() {
  DATA_OBJECT factlist;
  EnvGetFactList(clips_environment, &factlist, NULL);
  std::cout << "we have " << GetDOLength(factlist) << " facts!\n";

  DATA_OBJECT valid_two_body_decay_trees;
  EnvEval(clips_environment,
      "(do-for-all-facts ((?f ValidTwoBodyDecayTrees) (?g TwoBodyDecayTree)) (member$ ?g ?f:two_body_decay_trees) (fact-slot-value ?f two_body_decay_trees))",
      &valid_two_body_decay_trees);
  std::cout << "we have " << GetDOLength(valid_two_body_decay_trees)
      << " decay trees!\n";
  //EnvGetFactSlot(clips_environment, result_fact, "two_body_decay_trees",
  //    &valid_two_body_decay_trees);
  for (unsigned int tree_index = 1;
      tree_index <= GetDOLength(valid_two_body_decay_trees); ++tree_index) {
    DATA_OBJECT two_body_decay_tree;
    EnvGetFactSlot(clips_environment,
        GetMFValue(GetValue(valid_two_body_decay_trees), tree_index),
        "two_body_decays", &two_body_decay_tree);

    //std::cout << "this tree contains " << GetDOLength(two_body_decay_tree)
    //    << " two body decays!\n";

    for (unsigned int two_body_decay_index = 1;
        two_body_decay_index <= GetDOLength(two_body_decay_tree);
        ++two_body_decay_index) {

      /*DATA_OBJECT blub2;
       EnvFactSlotNames(clips_environment,
       GetMFValue(GetValue(two_body_decay_tree), two_body_decay_index),
       &blub2);
       for (unsigned int name_index = 1; name_index <= GetDOLength(blub2);
       ++name_index) {
       std::cout << ValueToString(GetMFValue(GetValue(blub2), name_index))
       << std::endl;
       }*/

      DATA_OBJECT mother;
      EnvGetFactSlot(clips_environment,
          GetMFValue(GetValue(two_body_decay_tree), two_body_decay_index),
          "mother", &mother);
      DATA_OBJECT daughter1;
      EnvGetFactSlot(clips_environment,
          GetMFValue(GetValue(two_body_decay_tree), two_body_decay_index),
          "daughter1", &daughter1);
      DATA_OBJECT daughter2;
      EnvGetFactSlot(clips_environment,
          GetMFValue(GetValue(two_body_decay_tree), two_body_decay_index),
          "daughter2", &daughter2);

      //DATA_OBJECT blub;
      //EnvFactSlotNames(clips_environment, &mother, &blub);
      //std::cout << "mother fact: " << GetDOLength(blub) << " slots"
      //    << std::endl;
      //DATA_OBJECT particle;
      //EnvGetFactSlot(clips_environment, GetValue(mother_fact), "unique_id",
      //    &particle);
      //std::cout << "test me: " << *((char*) particle.value) << std::endl;

      ParticleStateInfo mother_particle = getSpinWaveFromClipsEnvironment(
          ValueToInteger(GetValue(mother)));
      print(mother_particle);
      std::cout << "decays to\n";
      ParticleStateInfo daughter1_particle = getSpinWaveFromClipsEnvironment(
          ValueToInteger(GetValue(daughter1)));
      print(daughter1_particle);
      std::cout << " and ";
      ParticleStateInfo daughter2_particle = getSpinWaveFromClipsEnvironment(
          ValueToInteger(GetValue(daughter2)));
      print(daughter2_particle);
    }
  }
}

ParticleStateInfo DecayGenerator::getSpinWaveFromClipsEnvironment(
    int unique_id) const {
  ParticleStateInfo particle;

  DATA_OBJECT spin_waves;
  std::stringstream clips_query;
  clips_query << "(find-fact ((?m SpinWave)) (= ?m:unique_id " << unique_id
      << "))";

  EnvEval(clips_environment, clips_query.str().c_str(), &spin_waves);
  void *spin_wave_fact = GetMFValue(GetValue(spin_waves), 1);

  DATA_OBJECT data;

  // get particle id or id quantum numbers
  EnvGetFactSlot(clips_environment, spin_wave_fact, "unique_id", &data);
  particle.unique_id_ = ValueToInteger(GetValue(data));

  setPIDInfo(particle, spin_wave_fact);

  // get spin information
  EnvGetFactSlot(clips_environment, spin_wave_fact, "spin_num", &data);
  particle.spin_information_.J_numerator_ = ValueToInteger(GetValue(data));
  EnvGetFactSlot(clips_environment, spin_wave_fact, "spin_denom", &data);
  particle.spin_information_.J_denominator_ = ValueToInteger(GetValue(data));
  EnvGetFactSlot(clips_environment, spin_wave_fact, "spin_z_num", &data);
  particle.spin_information_.J_z_numerator_ = ValueToInteger(GetValue(data));
  return particle;
}

void DecayGenerator::setPIDInfo(ParticleStateInfo& particle,
    void* spin_wave_fact) const {
  // search for a particle in the list, that has the same QNs a the spin wave
  ComPWA::ParticleProperties particle_properties;
  DATA_OBJECT data;
  EnvGetFactSlot(clips_environment, spin_wave_fact, "charge", &data);
  particle_properties.charge_ = ValueToInteger(GetValue(data));

  EnvGetFactSlot(clips_environment, spin_wave_fact, "isospin_num", &data);
  particle_properties.isospin_.J_numerator_ = ValueToInteger(GetValue(data));
  EnvGetFactSlot(clips_environment, spin_wave_fact, "isospin_denom", &data);
  particle_properties.isospin_.J_denominator_ = ValueToInteger(GetValue(data));
  EnvGetFactSlot(clips_environment, spin_wave_fact, "isospin_z_num", &data);
  particle_properties.isospin_.J_z_numerator_ = ValueToInteger(GetValue(data));

  EnvGetFactSlot(clips_environment, spin_wave_fact, "spin_num", &data);
  particle_properties.spin_.J_numerator_ = ValueToInteger(GetValue(data));
  EnvGetFactSlot(clips_environment, spin_wave_fact, "spin_denom", &data);
  particle_properties.spin_.J_denominator_ = ValueToInteger(GetValue(data));
  EnvGetFactSlot(clips_environment, spin_wave_fact, "spin_z_num", &data);
  particle_properties.spin_.J_z_numerator_ = ValueToInteger(GetValue(data));

  EnvGetFactSlot(clips_environment, spin_wave_fact, "parity", &data);
  particle_properties.parity_ = ValueToInteger(GetValue(data));

  EnvGetFactSlot(clips_environment, spin_wave_fact, "cparity", &data);
  particle_properties.cparity_ = ValueToInteger(GetValue(data));

  std::vector<ComPWA::ParticleProperties> candidate_list =
      ComPWA::PhysConst::Instance().findParticlesWithQN(particle_properties);

  std::cout << "found " << candidate_list.size() << " candidates\n";
  for (unsigned int i = 0; i < candidate_list.size(); ++i) {
    std::cout << candidate_list[i].name_ << std::endl;
  }
}

void DecayGenerator::print(const ParticleStateInfo &psi) const {
  std::cout << "id: " << psi.unique_id_ << "; J="
      << psi.spin_information_.J_numerator_ << "/"
      << psi.spin_information_.J_denominator_ << "; J_z="
      << psi.spin_information_.J_z_numerator_ << "/"
      << psi.spin_information_.J_denominator_ << std::endl;
}



void DecayGenerator::createAllValidTwoBodyDecayTrees() {

}

/*std::vector<ParticleIndexDecayTree> DecayGenerator::makeConcreteDecayTrees(
 const std::vector<ParticleIndexDecayTree>& decay_tree_realization_templates) {
 std::vector<ParticleIndexDecayTree> concrete_decay_trees;

 return concrete_decay_trees;
 }*/

boost::property_tree::ptree DecayGenerator::portPropertyTree(
    const std::vector<ParticleIndexDecayTree>& concrete_decay_trees) const {
  return boost::property_tree::ptree();
}

} /* namespace HelicityFormalism */
