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

#include <chrono>

#include "Core/PhysConst.hpp"

#include "Physics/DecayTree/DecayGenerator.hpp"

// this include has to be after the compwa includes
// because it defines an INTEGER in the global scope that introduces clashes
#include "clips.h"

namespace ComPWA {
namespace Physics {
namespace DecayTree {

DecayGenerator::DecayGenerator() :
    unique_spin_state_counter_(1) {
  clips_environment_ = CreateEnvironment();

  allowed_particle_names_.push_back("gamma");
  allowed_particle_names_.push_back("pi0");
  allowed_particle_names_.push_back("f0_980");
  allowed_particle_names_.push_back("f0_1370");
  allowed_particle_names_.push_back("f2_1270");
  allowed_particle_names_.push_back("omega");
  allowed_particle_names_.push_back("jpsi");
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

    bool is_massless(false);
    if (0.0
        == ComPWA::PhysConst::Instance().findParticle(
            particle_info.particle_info_.name_).mass_)
      is_massless = true;

    for (int i = -1 * s.J_numerator_; i <= (int) s.J_numerator_;
        i = i + s.J_denominator_) {
      if (!(is_massless && i == 0)) {
        s.J_z_numerator_ = i;
        particle_info.spin_z_components_.push_back(s);
      }
    }
  }
}

int DecayGenerator::lookupPID(const std::string& name) const {
  return ComPWA::PhysConst::Instance().findParticle(name).id_;
}

DecayConfiguration DecayGenerator::createDecayConfiguration() {
  std::vector<SpinWaveDecayTree> two_body_decay_trees = generate();

  DecayConfiguration decay_configuration;

  initializeAllowedParticlePool();

  for (unsigned int i = 0; i < two_body_decay_trees.size(); ++i) {
    addSpinWaveTwoBodyDecayToDecayConfiguration(decay_configuration,
        two_body_decay_trees[i]);
  }

  return decay_configuration;
}

void DecayGenerator::initializeAllowedParticlePool() {
  for (unsigned int i = 0; i < allowed_particle_names_.size(); ++i) {
    allowed_particle_pool_.push_back(
        ComPWA::PhysConst::Instance().findParticle(allowed_particle_names_[i]));
  }
}

void DecayGenerator::addSpinWaveTwoBodyDecayToDecayConfiguration(
    DecayConfiguration& decay_configuration,
    const SpinWaveDecayTree& two_body_decay_tree) const {

  std::map<unsigned int, std::vector<ParticleStateInfo> > temp_state_pool;
  std::vector<
      std::vector<std::pair<ParticleStateInfo, std::vector<ParticleStateInfo> > > > decay_trees;

  for (auto const &decay_node : two_body_decay_tree.unique_decay_node_index_tree_) {
    if (temp_state_pool.find(decay_node.first) == temp_state_pool.end()) {
      temp_state_pool[decay_node.first] =
          createParticleStateInfoCandidates(
              two_body_decay_tree.unique_decay_node_index_to_spin_wave_index_mapping_.at(
                  decay_node.first));
    }

    std::vector<std::vector<ParticleStateInfo> > daughter_list_combinations;
    for (unsigned int i = 0; i < decay_node.second.size(); ++i) {
      if (temp_state_pool.find(decay_node.second[i]) == temp_state_pool.end()) {
        temp_state_pool[decay_node.second[i]] =
            createParticleStateInfoCandidates(
                two_body_decay_tree.unique_decay_node_index_to_spin_wave_index_mapping_.at(
                    decay_node.second[i]));
      }

      if (daughter_list_combinations.size() == 0) {
        for (auto particle : temp_state_pool[decay_node.second[i]]) {
          std::vector<ParticleStateInfo> blank_vector;
          blank_vector.push_back(particle);
          daughter_list_combinations.push_back(blank_vector);
        }
      }
      else {
        std::vector<std::vector<ParticleStateInfo> > new_daughter_list_combinations;
        for (auto const& combination : daughter_list_combinations) {
          for (auto const& particle : temp_state_pool[decay_node.second[i]]) {
            std::vector<ParticleStateInfo> extended_combination(combination);
            extended_combination.push_back(particle);
            new_daughter_list_combinations.push_back(extended_combination);
          }
        }
        daughter_list_combinations = new_daughter_list_combinations;
      }
    }

    // add them to the decay trees
    for (auto const& mother : temp_state_pool[decay_node.first]) {
      for (auto const& daughters : daughter_list_combinations) {
        if (decay_trees.size() == 0) {
          std::vector<
              std::pair<ParticleStateInfo, std::vector<ParticleStateInfo> > > blank_decay_tree;
          auto decay_pair = std::make_pair(mother, daughters);
          blank_decay_tree.push_back(decay_pair);
          decay_trees.push_back(blank_decay_tree);
        }
        else {
          std::vector<
              std::vector<
                  std::pair<ParticleStateInfo, std::vector<ParticleStateInfo> > > > new_decay_trees;
          for (auto const& current_decay_tree : decay_trees) {
            std::vector<
                std::pair<ParticleStateInfo, std::vector<ParticleStateInfo> > > extended_decay_tree(
                current_decay_tree);
            auto decay_pair = std::make_pair(mother, daughters);
            extended_decay_tree.push_back(decay_pair);
            new_decay_trees.push_back(extended_decay_tree);
          }
          decay_trees = new_decay_trees;
        }
      }
    }
  }

  for (auto const& decay_tree : decay_trees) {
    // for each possible mother state make copy for all pairs of daughters
    for (auto const &decay_node : decay_tree) {
      const boost::property_tree::ptree decay_strength_info_and_phase;

      decay_configuration.addDecayToCurrentDecayTree(decay_node.first,
          decay_node.second, decay_strength_info_and_phase);
    }
    decay_configuration.addCurrentDecayTreeToList();
  }
}

std::vector<ParticleStateInfo> DecayGenerator::createParticleStateInfoCandidates(
    unsigned int spin_wave_index) const {
  SpinWave spin_wave = all_spin_waves_[spin_wave_index];

  std::vector<ParticleStateInfo> candidates;

  // find candidates
  for (unsigned int i = 0; i < allowed_particle_pool_.size(); ++i) {
    Spin wave_spin =
        spin_wave.spin_like_quantum_numbers_[ComPWA::PhysConst::Instance().getQuantumNumberName(
            ComPWA::QuantumNumbers::SPIN)];
    if (wave_spin.J_numerator_ / wave_spin.J_denominator_
        == allowed_particle_pool_[i].spin_.J_numerator_
            / allowed_particle_pool_[i].spin_.J_denominator_) {
      if (wave_spin.J_z_numerator_ == 0
          && 0.0 == allowed_particle_pool_[i].mass_) {
        continue;
      }
      else {
        ParticleStateInfo ps;
        ps.spin_information_ = wave_spin;

        ps.pid_information_.particle_id_ = allowed_particle_pool_[i].id_;
        ps.pid_information_.name_ = allowed_particle_pool_[i].name_;

        ps.dynamical_information_ = createDynamicInfo(allowed_particle_pool_[i],
            DynamicalFunctions::DynamicalInfoTypes::RELATIVE_BREIT_WIGNER);

        candidates.push_back(ps);
      }
    }
  }

  return candidates;
}

DynamicalInfo DecayGenerator::createDynamicInfo(
    const ParticleProperties& particle_properties,
    ComPWA::Physics::DynamicalFunctions::DynamicalInfoTypes dynamical_type) const {
  DynamicalInfo dynamical_info;

  if (dynamical_type
      == ComPWA::Physics::DynamicalFunctions::DynamicalInfoTypes::RELATIVE_BREIT_WIGNER) {
    dynamical_info.put("DynamicalInfo.type",
        ComPWA::Physics::DynamicalFunctions::DynamicalTypeToString.at(
            dynamical_type));
    dynamical_info.put("DynamicalInfo.mass", particle_properties.mass_);
    dynamical_info.put("DynamicalInfo.mass_fix", 1);
    dynamical_info.put("DynamicalInfo.mass_min",
        0.5 * particle_properties.mass_);
    dynamical_info.put("DynamicalInfo.mass_max",
        1.5 * particle_properties.mass_);
    dynamical_info.put("DynamicalInfo.width", particle_properties.width_);
    dynamical_info.put("DynamicalInfo.width_fix", 1);
    dynamical_info.put("DynamicalInfo.width_min",
        0.5 * particle_properties.width_);
    dynamical_info.put("DynamicalInfo.width_max",
        10.0 * particle_properties.width_);
    dynamical_info.put("DynamicalInfo.mesonRadius", 1.0);
    dynamical_info.put("DynamicalInfo.norm", 1);
    dynamical_info.put("DynamicalInfo.par1", 1.0);
    dynamical_info.put("DynamicalInfo.par2", 1.0);
  }

  return dynamical_info;
}

std::vector<SpinWaveDecayTree> DecayGenerator::generate() {
  setupClipsEnvironment();

  /*createAllSpinWaves();
   createAllValidTwoBodyDecays();
   std::cout << "number of spin waves: " << all_spin_waves_.size() << std::endl;
   //std::cout << "number of two body decays: "
   //    << all_spin_wave_two_body_decays_.size() << std::endl;
   *
   *
   */

  addAllowedQuantumNumbersToClipsEnvironment();
  addConservedQuantumNumbersToClipsEnvironment();
  addAllInitialAndFinalStateCombinationsToClipsEnvironment();

// run
//EnvFacts(clips_environment_, "stdout", NULL, -1, -1, -1);
  int max_number_of_rules_to_fire(-1);
  EnvRun(clips_environment_, max_number_of_rules_to_fire);
//EnvFacts(clips_environment_, "stdout", NULL, -1, -1, -1);

// get expertise
  std::vector<SpinWaveDecayTree> two_body_decay_trees = getExpertise();

//high_resolution_clock::time_point t4 = high_resolution_clock::now();
  EnvReset(clips_environment_);

  std::cout << "we have " << two_body_decay_trees.size() << " decay trees!"
      << std::endl;

  for (unsigned int i = 0; i < two_body_decay_trees.size(); ++i) {
    printDecayTree(two_body_decay_trees[i]);
  }

  return two_body_decay_trees;
}

/*void DecayGenerator::createAllSpinWaves() {
 unsigned int unique_id_counter(0);
 SpinWave spin_wave;
 for (unsigned int spin_index = 0; spin_index < allowed_spins_.size();
 ++spin_index) {
 spin_wave.spin_ = allowed_spins_[spin_index];
 for (unsigned int isospin_index = 0;
 isospin_index < allowed_isospins_.size(); ++isospin_index) {
 spin_wave.isospin_ = allowed_isospins_[isospin_index];
 for (unsigned int charge_index = 0;
 charge_index < allowed_charges_.size(); ++charge_index) {
 spin_wave.charge_ = allowed_charges_[charge_index];
 for (unsigned int parity_index = 0;
 parity_index < allowed_parities_.size(); ++parity_index) {
 spin_wave.parity_ = allowed_parities_[parity_index];
 for (unsigned int cparity_index = 0;
 cparity_index < allowed_cparites_.size(); ++cparity_index) {
 spin_wave.cparity_ = allowed_cparites_[cparity_index];
 spin_wave.unique_id_ = ++unique_id_counter;
 all_spin_waves_.push_back(spin_wave);
 }
 }
 }
 }
 }
 }

 void DecayGenerator::createAllValidTwoBodyDecays() {
 // loop over all spin waves
 SpinWaveTwoBodyDecay swtbd;
 for (unsigned int mother_id = 0; mother_id < all_spin_waves_.size();
 ++mother_id) {
 swtbd.mother_index_ = mother_id;
 std::cout << mother_id << std::endl;
 for (unsigned int daughter1_id = 0; daughter1_id < all_spin_waves_.size();
 ++daughter1_id) {
 for (unsigned int daughter2_id = 0; daughter2_id < all_spin_waves_.size();
 ++daughter2_id) {
 swtbd.daughter_indices_ = std::make_pair(daughter1_id, daughter2_id);
 if (validateTwoBodyDecay(swtbd)) {
 IndexList daughter_indices;
 daughter_indices.push_back(daughter1_id);
 daughter_indices.push_back(daughter2_id);
 std::sort(daughter_indices.begin(), daughter_indices.end());
 all_spin_wave_two_body_decays_.push_back(swtbd);
 two_body_decays_lookup_table_[daughter_indices].push_back(
 all_spin_wave_two_body_decays_.size() - 1);
 }
 }
 }
 }
 }

 bool DecayGenerator::validateTwoBodyDecay(
 const SpinWaveTwoBodyDecay& two_body_decay) {
 //  EnvIncrementGCLocks(clips_environment);

 //using namespace std::chrono;

 //high_resolution_clock::time_point t1 = high_resolution_clock::now();

 addTwoBodyDecayToClipsEnviroment(two_body_decay);

 //high_resolution_clock::time_point t2 = high_resolution_clock::now();

 // run
 //EnvFacts(clips_environment, "stdout", NULL, -1, -1, -1);
 int max_number_of_rules_to_fire(-1);
 EnvRun(clips_environment_, max_number_of_rules_to_fire);
 //EnvFacts(clips_environment, "stdout", NULL, -1, -1, -1);

 //high_resolution_clock::time_point t3 = high_resolution_clock::now();

 // get expertise
 std::vector<int> violated_rules_ids = getExpertise();

 //high_resolution_clock::time_point t4 = high_resolution_clock::now();

 EnvReset(clips_environment_);

 //high_resolution_clock::time_point t5 = high_resolution_clock::now();

 /*duration<double> time_span1 = duration_cast<duration<double>>(t2 - t1);
 duration<double> time_span2 = duration_cast<duration<double>>(t3 - t2);
 duration<double> time_span3 = duration_cast<duration<double>>(t4 - t3);
 duration<double> time_span4 = duration_cast<duration<double>>(t5 - t4);

 std::cout<<"adding two body decay takes (seconds): "<< time_span1.count() <<std::endl;
 std::cout<<"running clips (seconds): "<< time_span2.count() <<std::endl;
 std::cout<<"getting results (seconds): "<< time_span3.count() <<std::endl;
 std::cout<<"resetting (seconds): "<< time_span4.count() <<std::endl;*/
/*
 bool result(true);

 if (violated_rules_ids.size() > 0) {
 //std::cout << "this decay violates " << violated_rules_ids.size()
 //    << " rules\n";
 result = false;
 }

 return result;
 }*/

void DecayGenerator::setupClipsEnvironment() const {
  std::stringstream path;

  // first load the general file that includes all the templates
  // and important runtime functions
  path << getenv("COMPWA_DIR") << "/Physics/DecayTree/General.clp";
  EnvLoad(clips_environment_, path.str().c_str());
  path.str("");

  // then load all the decay conservation laws for your decay
  path << getenv("COMPWA_DIR")
      << "/Physics/DecayTree/StrictlyConservedQNRules.clp";
  EnvLoad(clips_environment_, path.str().c_str());

  path.str("");
  path << getenv("COMPWA_DIR") << "/Physics/DecayTree/Parity.clp";
  EnvLoad(clips_environment_, path.str().c_str());

  path.str("");
  path << getenv("COMPWA_DIR") << "/Physics/DecayTree/HelicityFormalism.clp";
  EnvLoad(clips_environment_, path.str().c_str());

  path.str("");

  // this has always to be at the end
  path << getenv("COMPWA_DIR")
      << "/Physics/DecayTree/TwoBodyDecayTreeGenerator.clp";
  EnvLoad(clips_environment_, path.str().c_str());

  EnvReset(clips_environment_);
}

void DecayGenerator::addAllowedQuantumNumbersToClipsEnvironment() {
  for (auto const& allowed_spin_like_quantum_number : allowed_spin_like_quantum_numbers_) {
    addSpinLikeAllowedQuantumNumbersToClipsEnvironment(
        allowed_spin_like_quantum_number);
  }

  for (auto const& allowed_integer_like_quantum_number : allowed_integer_like_quantum_numbers_) {
    addIntLikeAllowedQuantumNumbersToClipsEnvironment(
        allowed_integer_like_quantum_number);
  }
}

void DecayGenerator::addSpinLikeAllowedQuantumNumbersToClipsEnvironment(
    const AllowedQuantumNumbers<ComPWA::Spin>& allowed_qn) {
  std::vector<int> unique_spin_state_list;

// loop over spin states (SpinQuantumNumber (unique_id 1) (numerator 0) (denominator 1) (z_component_numerator 0))
  for (auto spin_state = allowed_qn.allowed_values_.begin();
      spin_state != allowed_qn.allowed_values_.end(); ++spin_state) {
    unique_spin_state_list.push_back(
        addSpinQuantumNumberToClipsEnvironment(*spin_state));
  }

  AllowedQuantumNumbers<int> dummy;
  dummy.allowed_values_ = unique_spin_state_list;
  dummy.quantum_number_name_ = allowed_qn.quantum_number_name_;
  dummy.type = allowed_qn.type;
  dummy.required_quantum_numbers_names_ = allowed_qn.required_quantum_numbers_names_;
  addIntLikeAllowedQuantumNumbersToClipsEnvironment(dummy);
}

unsigned int DecayGenerator::addSpinQuantumNumberToClipsEnvironment(
    const ComPWA::Spin& spin_state) {
  int return_value(-1);
// check if this spin state exists otherwise create and add unique spin state index to list
  DATA_OBJECT found_spin_quantum_number_facts;
  std::stringstream clips_query;
  clips_query << "(find-fact ((?f SpinQuantumNumber)) (and (= ?f:numerator "
      << spin_state.J_numerator_ << ")"
          " (and (= ?f:denominator " << spin_state.J_denominator_ << ") "
          "(= ?f:z_component_numerator " << spin_state.J_z_numerator_ << "))))";

  /*clips_query << "(find-spin-quantum-number-list " << spin_state.J_numerator_
   << " " << spin_state.J_denominator_ << " " << spin_state.J_z_numerator_
   << ")";*/

  EnvEval(clips_environment_, clips_query.str().c_str(),
      &found_spin_quantum_number_facts);

  if (0 < GetDOLength(found_spin_quantum_number_facts)) {
    DATA_OBJECT data;
    EnvGetFactSlot(clips_environment_,
        GetMFValue(GetValue(found_spin_quantum_number_facts), 1), "unique_id",
        &data);
    return_value = ValueToInteger(GetValue(data));
  }
  else {
    void* spin_qn_template = EnvFindDeftemplate(clips_environment_,
        "SpinQuantumNumber");
    void* spin_qn_fact = EnvCreateFact(clips_environment_, spin_qn_template);
    if (spin_qn_fact != NULL) {
      DATA_OBJECT field;
      field.type = INTEGER;

      field.value = EnvAddLong(clips_environment_, unique_spin_state_counter_);
      EnvPutFactSlot(clips_environment_, spin_qn_fact, "unique_id", &field);
      field.value = EnvAddLong(clips_environment_, spin_state.J_numerator_);
      EnvPutFactSlot(clips_environment_, spin_qn_fact, "numerator", &field);
      field.value = EnvAddLong(clips_environment_, spin_state.J_denominator_);
      EnvPutFactSlot(clips_environment_, spin_qn_fact, "denominator", &field);
      field.value = EnvAddLong(clips_environment_, spin_state.J_z_numerator_);
      EnvPutFactSlot(clips_environment_, spin_qn_fact, "z_component_numerator",
          &field);

      EnvAssert(clips_environment_, spin_qn_fact);
      return_value = unique_spin_state_counter_;
      ++unique_spin_state_counter_;
    }
  }
  return return_value;
}

void DecayGenerator::addIntLikeAllowedQuantumNumbersToClipsEnvironment(
    const AllowedQuantumNumbers<int>& allowed_qn) const {
// add the allowed qn fact to clips
  std::string name;
  if (allowed_qn.type == QuantumNumberTypes::COMPOSITE_PARTICLE_BASED)
    name = "AllowedDecayQuantumNumbers";
  else
    name = "AllowedQuantumNumbers";
  void* allowed_spin_qn_template = EnvFindDeftemplate(clips_environment_,
      name.c_str());
  void* allowed_spin_qn_fact = EnvCreateFact(clips_environment_,
      allowed_spin_qn_template);
  if (allowed_spin_qn_fact != NULL) {
    DATA_OBJECT field;
    field.type = STRING;

    field.value = EnvAddSymbol(clips_environment_,
        allowed_qn.quantum_number_name_.c_str());
    EnvPutFactSlot(clips_environment_, allowed_spin_qn_fact, "name", &field);

    DATA_OBJECT allowed_values;
    void* multifield_ptr = EnvCreateMultifield(clips_environment_,
        allowed_qn.allowed_values_.size());

    for (unsigned int i = 0; i < allowed_qn.allowed_values_.size(); ++i) {
      SetMFType(multifield_ptr, i + 1, INTEGER);
      SetMFValue(multifield_ptr, i + 1,
          EnvAddLong(clips_environment_, allowed_qn.allowed_values_[i]));
    }
    SetType(allowed_values, MULTIFIELD);
    SetValue(allowed_values, multifield_ptr);

    SetDOBegin(allowed_values, 1);
    SetDOEnd(allowed_values, allowed_qn.allowed_values_.size());

    EnvPutFactSlot(clips_environment_, allowed_spin_qn_fact, "values",
        &allowed_values);

    DATA_OBJECT required_names;
    multifield_ptr = EnvCreateMultifield(clips_environment_,
        allowed_qn.required_quantum_numbers_names_.size());

    for (unsigned int i = 0;
        i < allowed_qn.required_quantum_numbers_names_.size(); ++i) {
      SetMFType(multifield_ptr, i + 1, STRING);
      SetMFValue(multifield_ptr, i + 1,
          EnvAddSymbol(clips_environment_,
              allowed_qn.required_quantum_numbers_names_[i].c_str()));
    }
    SetType(required_names, MULTIFIELD);
    SetValue(required_names, multifield_ptr);

    SetDOBegin(required_names, 1);
    SetDOEnd(required_names, allowed_qn.required_quantum_numbers_names_.size());

    EnvPutFactSlot(clips_environment_, allowed_spin_qn_fact,
        "required_variable_names", &required_names);

    EnvAssert(clips_environment_, allowed_spin_qn_fact);
  }
}

void DecayGenerator::addConservedQuantumNumbersToClipsEnvironment() {
  void* name_list_template = EnvFindDeftemplate(clips_environment_, "NameList");

  void* name_list_fact = EnvCreateFact(clips_environment_, name_list_template);
  if (name_list_fact != NULL) {
    DATA_OBJECT field;
    field.type = STRING;

    field.value = EnvAddSymbol(clips_environment_, "ConservedQuantumNumbers");
    EnvPutFactSlot(clips_environment_, name_list_fact, "name", &field);

    DATA_OBJECT qn_names;
    void* multifield_ptr = EnvCreateMultifield(clips_environment_,
        conserved_quantum_numbers_.size());

    for (unsigned int i = 0; i < conserved_quantum_numbers_.size(); ++i) {
      SetMFType(multifield_ptr, i + 1, STRING);
      SetMFValue(multifield_ptr, i + 1,
          EnvAddSymbol(clips_environment_,
              conserved_quantum_numbers_[i].c_str()));
    }
    SetType(qn_names, MULTIFIELD);
    SetValue(qn_names, multifield_ptr);

    SetDOBegin(qn_names, 1);
    SetDOEnd(qn_names, conserved_quantum_numbers_.size());

    EnvPutFactSlot(clips_environment_, name_list_fact, "names", &qn_names);

    EnvAssert(clips_environment_, name_list_fact);
  }
}

void DecayGenerator::addAllInitialAndFinalStateCombinationsToClipsEnvironment() {
  std::vector<std::pair<SpinWave, std::vector<SpinWave> > > IF_state_combinations;
  std::vector<std::pair<SpinWave, std::vector<SpinWave> > > temp_IF_state_combinations;

// create all initial and final state combinations first
  if (final_state_particles_.size() > 0) {
    //initial state combinations first

    for (unsigned int i = 0;
        i < mother_state_particle_.spin_z_components_.size(); ++i) {
      IF_state_combinations.push_back(
          std::make_pair(createSpinWave(mother_state_particle_, i),
              std::vector<SpinWave>()));
    }

    for (unsigned int fsp_index = 0; fsp_index < final_state_particles_.size();
        ++fsp_index) {
      for (unsigned int i = 0;
          i < final_state_particles_[fsp_index].spin_z_components_.size();
          ++i) {
        for (unsigned int if_combination_index = 0;
            if_combination_index < IF_state_combinations.size();
            ++if_combination_index) {
          std::pair<SpinWave, std::vector<SpinWave> > current_combination =
              IF_state_combinations[if_combination_index];
          current_combination.second.push_back(
              createSpinWave(final_state_particles_[fsp_index], i));
          temp_IF_state_combinations.push_back(current_combination);
        }
      }
      IF_state_combinations = temp_IF_state_combinations;
      temp_IF_state_combinations.clear();
    }

    // now add all the initial and final states to clips
    for (unsigned int i = 0; i < IF_state_combinations.size(); ++i) {
      addInitialAndFinalStateToClipsEnvironment(IF_state_combinations[i].first,
          IF_state_combinations[i].second);
    }
  }
}

SpinWave DecayGenerator::createSpinWave(const IFParticleInfo& particle,
    unsigned int state_index) const {
  SpinWave wave;

  ComPWA::ParticleProperties particle_properties =
      ComPWA::PhysConst::Instance().findParticle(particle.particle_info_.name_);

  wave.spin_like_quantum_numbers_[ComPWA::PhysConst::Instance().getQuantumNumberName(
      ComPWA::QuantumNumbers::SPIN)] = particle.spin_z_components_[state_index];
  wave.spin_like_quantum_numbers_[ComPWA::PhysConst::Instance().getQuantumNumberName(
      ComPWA::QuantumNumbers::ISOSPIN)] = particle_properties.isospin_;
  wave.integer_like_quantum_numbers_[ComPWA::PhysConst::Instance().getQuantumNumberName(
      ComPWA::QuantumNumbers::CHARGE)] = particle_properties.charge_;
  wave.integer_like_quantum_numbers_[ComPWA::PhysConst::Instance().getQuantumNumberName(
      ComPWA::QuantumNumbers::PARITY)] = particle_properties.parity_;
  wave.integer_like_quantum_numbers_[ComPWA::PhysConst::Instance().getQuantumNumberName(
      ComPWA::QuantumNumbers::CPARITY)] = particle_properties.cparity_;

  return wave;
}

void DecayGenerator::addInitialAndFinalStateToClipsEnvironment(
    const SpinWave& initial_state, const std::vector<SpinWave>& final_state) {
  void* ifs_template = EnvFindDeftemplate(clips_environment_,
      "InitialAndFinalState");
  EnvIncrementGCLocks(clips_environment_);
  void* ifs_fact = EnvCreateFact(clips_environment_, ifs_template);
  if (ifs_fact != NULL) {
    DATA_OBJECT field;
    field.type = FACT_ADDRESS;

    field.value =    //EnvAssert(clips_environment_,
        addSpinWaveToClipsEnvironment(initial_state);
    EnvPutFactSlot(clips_environment_, ifs_fact, "initial_state", &field);


    void* multifield_ptr = EnvCreateMultifield(clips_environment_,
        final_state.size());

    for (unsigned int i = 0; i < final_state.size(); ++i) {
      SetMFType(multifield_ptr, i + 1, FACT_ADDRESS);
      SetMFValue(multifield_ptr, i + 1,
          addSpinWaveToClipsEnvironment(final_state[i]));
    }

    DATA_OBJECT final_states;
    SetDOBegin(final_states, 1);
    SetDOEnd(final_states, final_state.size());

    SetType(final_states, MULTIFIELD);
    SetValue(final_states, multifield_ptr);

    EnvPutFactSlot(clips_environment_, ifs_fact, "final_state", &final_states);
    EnvAssignFactSlotDefaults(clips_environment_, ifs_fact);
    EnvAssert(clips_environment_, ifs_fact);
  }
  EnvDecrementGCLocks(clips_environment_);
}

void* DecayGenerator::addSpinWaveToClipsEnvironment(const SpinWave& spinwave) {
  std::vector<std::pair<std::string, int> > spin_qn_name_value_pairs;
  for (auto spin_like_qn = spinwave.spin_like_quantum_numbers_.begin();
      spin_like_qn != spinwave.spin_like_quantum_numbers_.end();
      ++spin_like_qn) {
    spin_qn_name_value_pairs.push_back(
        std::make_pair(spin_like_qn->first,
            addSpinQuantumNumberToClipsEnvironment(spin_like_qn->second)));
  }

  std::stringstream clips_query;
  std::stringstream values_part;
  values_part << "(explode$ \"";
  clips_query << "(find-spinwave-fact-list (explode$ \"";
  for (auto spin_like_qn = spin_qn_name_value_pairs.begin();
      spin_like_qn != spin_qn_name_value_pairs.end(); ++spin_like_qn) {
    clips_query << "\\\"" << spin_like_qn->first << "\\\" ";
    values_part << spin_like_qn->second << " ";
  }
  for (auto int_like_qn = spinwave.integer_like_quantum_numbers_.begin();
      int_like_qn != spinwave.integer_like_quantum_numbers_.end();
      ++int_like_qn) {
    clips_query << "\\\"" << int_like_qn->first << "\\\" ";
    values_part << int_like_qn->second << " ";
  }
  for (auto double_like_qn = spinwave.double_like_quantum_numbers_.begin();
      double_like_qn != spinwave.double_like_quantum_numbers_.end();
      ++double_like_qn) {
    clips_query << "\\\"" << double_like_qn->first << "\\\" ";
    values_part << double_like_qn->second << " ";
  }
  values_part << "\")";
  clips_query << "\") " << values_part.str() << ")";

  DATA_OBJECT found_spin_waves_facts;
  EnvEval(clips_environment_, clips_query.str().c_str(),
      &found_spin_waves_facts);

  void* spinwave_fact;
  if (0 < GetDOLength(found_spin_waves_facts)) {
    spinwave_fact = GetMFValue(GetValue(found_spin_waves_facts), 1);
  }
  else {
    void* spinwave_template = EnvFindDeftemplate(clips_environment_,
        "SpinWave");
    // set the facts
    spinwave_fact = EnvCreateFact(clips_environment_, spinwave_template);
    if (spinwave_fact != NULL) {

      unsigned int total_qn_count = spinwave.spin_like_quantum_numbers_.size()
          + spinwave.integer_like_quantum_numbers_.size()
          + spinwave.double_like_quantum_numbers_.size();

      void* qn_names_ptr = EnvCreateMultifield(clips_environment_,
          total_qn_count);

      void* qn_values_ptr = EnvCreateMultifield(clips_environment_,
          total_qn_count);

      unsigned int counter(1);
      for (auto spin_like_qn = spin_qn_name_value_pairs.begin();
          spin_like_qn != spin_qn_name_value_pairs.end(); ++spin_like_qn) {
        SetMFType(qn_names_ptr, counter, STRING);
        SetMFValue(qn_names_ptr, counter,
            EnvAddSymbol(clips_environment_, spin_like_qn->first.c_str()));
        SetMFType(qn_values_ptr, counter, INTEGER);
        SetMFValue(qn_values_ptr, counter,
            EnvAddLong(clips_environment_, spin_like_qn->second));
        ++counter;
      }
      for (auto int_like_qn = spinwave.integer_like_quantum_numbers_.begin();
          int_like_qn != spinwave.integer_like_quantum_numbers_.end();
          ++int_like_qn) {
        SetMFType(qn_names_ptr, counter, STRING);
        SetMFValue(qn_names_ptr, counter,
            EnvAddSymbol(clips_environment_, int_like_qn->first.c_str()));
        SetMFType(qn_values_ptr, counter, INTEGER);
        SetMFValue(qn_values_ptr, counter,
            EnvAddLong(clips_environment_, int_like_qn->second));
        ++counter;
      }
      for (auto double_like_qn = spinwave.double_like_quantum_numbers_.begin();
          double_like_qn != spinwave.double_like_quantum_numbers_.end();
          ++double_like_qn) {
        SetMFType(qn_names_ptr, counter, STRING);
        SetMFValue(qn_names_ptr, counter,
            EnvAddSymbol(clips_environment_, double_like_qn->first.c_str()));
        SetMFType(qn_values_ptr, counter, INTEGER);
        SetMFValue(qn_values_ptr, counter,
            EnvAddLong(clips_environment_, double_like_qn->second));
        ++counter;
      }

      DATA_OBJECT qn_names;
      DATA_OBJECT qn_values;

      SetType(qn_names, MULTIFIELD);
      SetValue(qn_names, qn_names_ptr);

      SetDOBegin(qn_names, 1);
      SetDOEnd(qn_names, total_qn_count);

      SetType(qn_values, MULTIFIELD);
      SetValue(qn_values, qn_values_ptr);

      SetDOBegin(qn_values, 1);
      SetDOEnd(qn_values, total_qn_count);

      EnvPutFactSlot(clips_environment_, spinwave_fact, "quantum_number_names",
          &qn_names);
      EnvPutFactSlot(clips_environment_, spinwave_fact, "quantum_number_values",
          &qn_values);

      EnvAssignFactSlotDefaults(clips_environment_, spinwave_fact);
      EnvAssert(clips_environment_, spinwave_fact);
    }
  }
  return spinwave_fact;
}

std::vector<SpinWaveDecayTree> DecayGenerator::getExpertise() {
  std::vector<SpinWaveDecayTree> valid_decay_trees;

  DATA_OBJECT decay_tree_fact_list;
  std::stringstream clips_query;
  clips_query << "(find-all-facts ((?f DecayTree)) TRUE)";
  EnvEval(clips_environment_, clips_query.str().c_str(), &decay_tree_fact_list);

  std::vector<void*> temp_fact_storage;

  for (unsigned int decay_tree_fact_index = GetDOBegin(decay_tree_fact_list);
      decay_tree_fact_index <= GetDOEnd(decay_tree_fact_list);
      ++decay_tree_fact_index) {
    temp_fact_storage.push_back(
        GetMFValue(GetValue(decay_tree_fact_list), decay_tree_fact_index));
  }
  for (unsigned int i = 0; i < temp_fact_storage.size(); ++i) {
    SpinWaveDecayTree swdt(
        getDecayTreeFromClipsEnvironment(temp_fact_storage[i]));
    if (std::find(valid_decay_trees.begin(), valid_decay_trees.end(), swdt)
        == valid_decay_trees.end()) {
      valid_decay_trees.push_back(swdt);
    }
  }

  return valid_decay_trees;
}

SpinWaveDecayTree DecayGenerator::getDecayTreeFromClipsEnvironment(
    void* decay_tree_fact) {
  SpinWaveDecayTree spin_wave_decay_tree;

  IndexList temp_indices;
  DATA_OBJECT all_occuring_wave_fact_list;
  EnvGetFactSlot(clips_environment_, decay_tree_fact, "all_occuring_waves",
      &all_occuring_wave_fact_list);
  for (unsigned int i = GetDOBegin(all_occuring_wave_fact_list);
      i <= GetDOEnd(all_occuring_wave_fact_list); ++i) {
    temp_indices.push_back(
        createSpinWave(GetMFValue(GetValue(all_occuring_wave_fact_list), i)));
  }

  DATA_OBJECT initial_state_wave_fact;
  EnvGetFactSlot(clips_environment_, decay_tree_fact, "initial_state_wave",
      &initial_state_wave_fact);
  spin_wave_decay_tree.top_node_unique_decay_node_index_ = createSpinWave(
      GetValue(initial_state_wave_fact));

  DATA_OBJECT unique_index_wave_mapping_fact_list;
  EnvGetFactSlot(clips_environment_, decay_tree_fact,
      "unique_index_wave_mapping", &unique_index_wave_mapping_fact_list);
  for (unsigned int i = GetDOBegin(unique_index_wave_mapping_fact_list);
      i <= GetDOEnd(unique_index_wave_mapping_fact_list); ++i) {
    DATA_OBJECT data;
    EnvGetFactSlot(clips_environment_,
        GetMFValue(GetValue(unique_index_wave_mapping_fact_list), i),
        "unique_id", &data);
    unsigned int unique_id = ValueToInteger(GetValue(data));
    EnvGetFactSlot(clips_environment_,
        GetMFValue(GetValue(unique_index_wave_mapping_fact_list), i),
        "list_index", &data);
    unsigned int list_index = ValueToInteger(GetValue(data));
    spin_wave_decay_tree.unique_decay_node_index_to_spin_wave_index_mapping_[unique_id] =
        temp_indices[list_index - 1];
  }

  DATA_OBJECT decay_fact_list;
  EnvGetFactSlot(clips_environment_, decay_tree_fact, "decays",
      &decay_fact_list);
  for (unsigned int i = GetDOBegin(decay_fact_list);
      i <= GetDOEnd(decay_fact_list); ++i) {
    DATA_OBJECT data;
    EnvGetFactSlot(clips_environment_, GetMFValue(GetValue(decay_fact_list), i),
        "mother", &data);
    unsigned int unique_id_mother = ValueToInteger(GetValue(data));

    IndexList daughter_ids;
    EnvGetFactSlot(clips_environment_, GetMFValue(GetValue(decay_fact_list), i),
        "daughters", &data);
    for (unsigned int daughter_index = GetDOBegin(data);
        daughter_index <= GetDOEnd(data); ++daughter_index) {
      daughter_ids.push_back(
          ValueToInteger(GetMFValue(GetValue(data), daughter_index)));
    }
    spin_wave_decay_tree.unique_decay_node_index_tree_[unique_id_mother] =
        daughter_ids;
  }
  return spin_wave_decay_tree;
}

unsigned int DecayGenerator::createSpinWave(void* spin_wave_fact) {
  SpinWave sw;

  DATA_OBJECT quantum_number_names_list;
  EnvGetFactSlot(clips_environment_, spin_wave_fact, "quantum_number_names",
      &quantum_number_names_list);

  DATA_OBJECT quantum_number_values_list;
  EnvGetFactSlot(clips_environment_, spin_wave_fact, "quantum_number_values",
      &quantum_number_values_list);

  for (unsigned int i = GetDOBegin(quantum_number_names_list);
      i <= GetDOEnd(quantum_number_names_list); ++i) {
    std::string name = ValueToString(
        GetMFValue(GetValue(quantum_number_names_list), i));
    if (name.find("spin") != std::string::npos) {
      sw.spin_like_quantum_numbers_[name] = createSpin(
          ValueToInteger(GetMFValue(GetValue(quantum_number_values_list), i)));
    }
    else {
      sw.integer_like_quantum_numbers_[name] = ValueToInteger(
          GetMFValue(GetValue(quantum_number_values_list), i));
    }
  }

  unsigned int position_index(all_spin_waves_.size());
  auto result = std::find_if(all_spin_waves_.begin(), all_spin_waves_.end(),
      sw);
  if (result == all_spin_waves_.end())
    all_spin_waves_.push_back(sw);
  else
    position_index = result - all_spin_waves_.begin();

//std::cout << sw << std::endl;
  return position_index;
}

ComPWA::Spin DecayGenerator::createSpin(
    int spin_quantum_number_unique_id) const {
  ComPWA::Spin spin;

  DATA_OBJECT spin_quantum_number_fact_list;
  std::stringstream clips_query;
  clips_query << "(find-fact ((?f SpinQuantumNumber)) (= ?f:unique_id "
      << spin_quantum_number_unique_id << "))";
  EnvEval(clips_environment_, clips_query.str().c_str(),
      &spin_quantum_number_fact_list);

//EnvFacts(clips_environment_, "stdout", NULL, -1, -1, -1);

  if (0 < GetDOLength(spin_quantum_number_fact_list)) {
    void* spin_quantum_number_fact = GetMFValue(
        GetValue(spin_quantum_number_fact_list),
        GetDOBegin(spin_quantum_number_fact_list));
    DATA_OBJECT data;
    EnvGetFactSlot(clips_environment_, spin_quantum_number_fact, "numerator",
        &data);
    spin.J_numerator_ = ValueToInteger(GetValue(data));
    EnvGetFactSlot(clips_environment_, spin_quantum_number_fact, "denominator",
        &data);
    spin.J_denominator_ = ValueToInteger(GetValue(data));
    EnvGetFactSlot(clips_environment_, spin_quantum_number_fact,
        "z_component_numerator", &data);
    spin.J_z_numerator_ = ValueToInteger(GetValue(data));
  }
  else {
    std::runtime_error(
        "DecayGenerator::createSpin: Error retrieving the spin fact from CLIPS, even though it is supposed to be there!");
  }
  return spin;
}

/*void DecayGenerator::addTwoBodyDecayToClipsEnviroment(
 const SpinWaveTwoBodyDecay& two_body_decay) const {

 // add the SpinWaves
 addSpinWaveToClipsEnvironment(all_spin_waves_[two_body_decay.mother_index_]);
 addSpinWaveToClipsEnvironment(
 all_spin_waves_[two_body_decay.daughter_indices_.first]);
 addSpinWaveToClipsEnvironment(
 all_spin_waves_[two_body_decay.daughter_indices_.second]);

 void* two_body_decay_template;
 two_body_decay_template = EnvFindDeftemplate(clips_environment_,
 "TwoBodyDecay");
 void* two_body_decay_fact = EnvCreateFact(clips_environment_,
 two_body_decay_template);
 if (two_body_decay_fact != NULL) {
 DATA_OBJECT field;
 field.type = INTEGER;

 field.value = EnvAddLong(clips_environment_,
 all_spin_waves_[two_body_decay.mother_index_].unique_id_);
 EnvPutFactSlot(clips_environment_, two_body_decay_fact, "mother", &field);
 field.value = EnvAddLong(clips_environment_,
 all_spin_waves_[two_body_decay.daughter_indices_.first].unique_id_);
 EnvPutFactSlot(clips_environment_, two_body_decay_fact, "daughter1",
 &field);
 field.value = EnvAddLong(clips_environment_,
 all_spin_waves_[two_body_decay.daughter_indices_.second].unique_id_);
 EnvPutFactSlot(clips_environment_, two_body_decay_fact, "daughter2",
 &field);

 EnvAssert(clips_environment_, two_body_decay_fact);
 }
 }

 void DecayGenerator::addSpinWaveToClipsEnvironment(
 const SpinWave& spin_wave) const {
 /*ComPWA::ParticleProperties particle_properties =
 ComPWA::PhysConst::Instance().findParticle(
 particle_info.particle_info_.name_);
 */
/*  void* spinwave_template;
 spinwave_template = EnvFindDeftemplate(clips_environment_, "SpinWave");
 // set the facts
 void* spinwave_fact = EnvCreateFact(clips_environment_, spinwave_template);
 if (spinwave_fact != NULL) {
 DATA_OBJECT data_value;
 data_value.type = INTEGER;
 data_value.value = EnvAddLong(clips_environment_, spin_wave.unique_id_);
 EnvPutFactSlot(clips_environment_, spinwave_fact, "unique_id", &data_value);
 data_value.value = EnvAddLong(clips_environment_, spin_wave.charge_);
 EnvPutFactSlot(clips_environment_, spinwave_fact, "charge", &data_value);
 data_value.value = EnvAddLong(clips_environment_,
 spin_wave.isospin_.J_numerator_);
 EnvPutFactSlot(clips_environment_, spinwave_fact, "isospin_num",
 &data_value);
 data_value.value = EnvAddLong(clips_environment_,
 spin_wave.isospin_.J_denominator_);
 EnvPutFactSlot(clips_environment_, spinwave_fact, "isospin_denom",
 &data_value);
 data_value.value = EnvAddLong(clips_environment_,
 spin_wave.isospin_.J_z_numerator_);
 EnvPutFactSlot(clips_environment_, spinwave_fact, "isospin_z_num",
 &data_value);
 data_value.value = EnvAddLong(clips_environment_,
 spin_wave.spin_.J_numerator_);
 EnvPutFactSlot(clips_environment_, spinwave_fact, "spin_num", &data_value);
 data_value.value = EnvAddLong(clips_environment_,
 spin_wave.spin_.J_denominator_);
 EnvPutFactSlot(clips_environment_, spinwave_fact, "spin_denom",
 &data_value);
 data_value.value = EnvAddLong(clips_environment_,
 spin_wave.spin_.J_z_numerator_);
 EnvPutFactSlot(clips_environment_, spinwave_fact, "spin_z_num",
 &data_value);

 data_value.value = EnvAddLong(clips_environment_, spin_wave.parity_);
 EnvPutFactSlot(clips_environment_, spinwave_fact, "parity", &data_value);
 data_value.value = EnvAddLong(clips_environment_, spin_wave.cparity_);
 EnvPutFactSlot(clips_environment_, spinwave_fact, "cparity", &data_value);

 EnvAssignFactSlotDefaults(clips_environment_, spinwave_fact);
 EnvAssert(clips_environment_, spinwave_fact);
 }
 }

 std::vector<int> DecayGenerator::getExpertise() {
 std::vector<int> violated_rules_id_codes;

 //DATA_OBJECT factlist;
 //EnvGetFactList(clips_environment, &factlist, NULL);
 //std::cout << "we have " << GetDOLength(factlist) << " facts!\n";

 DATA_OBJECT violated_rules_for_decay;
 std::stringstream clips_query;
 clips_query
 << "(do-for-all-facts ((?f ViolatingRulesForDecay)) (< -1 (length ?f:list_of_violated_rules)) (fact-slot-value ?f list_of_violated_rules))", EnvEval(
 clips_environment_, clips_query.str().c_str(), &violated_rules_for_decay);

 for (unsigned int violated_rule_index = 1;
 violated_rule_index <= GetDOLength(violated_rules_for_decay);
 ++violated_rule_index) {
 violated_rules_id_codes.push_back(
 ValueToInteger(
 GetMFValue(GetValue(violated_rules_for_decay), violated_rule_index)));
 }

 return violated_rules_id_codes;
 }

 ParticleStateInfo DecayGenerator::getSpinWaveFromClipsEnvironment(
 int unique_id) const {
 ParticleStateInfo particle;

 DATA_OBJECT spin_waves;
 std::stringstream clips_query;
 clips_query << "(find-fact ((?m SpinWave)) (= ?m:unique_id " << unique_id
 << "))";

 EnvEval(clips_environment_, clips_query.str().c_str(), &spin_waves);
 void *spin_wave_fact = GetMFValue(GetValue(spin_waves), 1);

 DATA_OBJECT data;

 // get particle id or id quantum numbers
 EnvGetFactSlot(clips_environment_, spin_wave_fact, "unique_id", &data);
 particle.unique_id_ = ValueToInteger(GetValue(data));

 setPIDInfo(particle, spin_wave_fact);

 // get spin information
 EnvGetFactSlot(clips_environment_, spin_wave_fact, "spin_num", &data);
 particle.spin_information_.J_numerator_ = ValueToInteger(GetValue(data));
 EnvGetFactSlot(clips_environment_, spin_wave_fact, "spin_denom", &data);
 particle.spin_information_.J_denominator_ = ValueToInteger(GetValue(data));
 EnvGetFactSlot(clips_environment_, spin_wave_fact, "spin_z_num", &data);
 particle.spin_information_.J_z_numerator_ = ValueToInteger(GetValue(data));
 return particle;
 }

 void DecayGenerator::setPIDInfo(ParticleStateInfo& particle,
 void* spin_wave_fact) const {
 // search for a particle in the list, that has the same QNs a the spin wave
 ComPWA::ParticleProperties particle_properties;
 DATA_OBJECT data;
 EnvGetFactSlot(clips_environment_, spin_wave_fact, "charge", &data);
 particle_properties.charge_ = ValueToInteger(GetValue(data));

 EnvGetFactSlot(clips_environment_, spin_wave_fact, "isospin_num", &data);
 particle_properties.isospin_.J_numerator_ = ValueToInteger(GetValue(data));
 EnvGetFactSlot(clips_environment_, spin_wave_fact, "isospin_denom", &data);
 particle_properties.isospin_.J_denominator_ = ValueToInteger(GetValue(data));
 EnvGetFactSlot(clips_environment_, spin_wave_fact, "isospin_z_num", &data);
 particle_properties.isospin_.J_z_numerator_ = ValueToInteger(GetValue(data));

 EnvGetFactSlot(clips_environment_, spin_wave_fact, "spin_num", &data);
 particle_properties.spin_.J_numerator_ = ValueToInteger(GetValue(data));
 EnvGetFactSlot(clips_environment_, spin_wave_fact, "spin_denom", &data);
 particle_properties.spin_.J_denominator_ = ValueToInteger(GetValue(data));
 EnvGetFactSlot(clips_environment_, spin_wave_fact, "spin_z_num", &data);
 particle_properties.spin_.J_z_numerator_ = ValueToInteger(GetValue(data));

 EnvGetFactSlot(clips_environment_, spin_wave_fact, "parity", &data);
 particle_properties.parity_ = ValueToInteger(GetValue(data));

 EnvGetFactSlot(clips_environment_, spin_wave_fact, "cparity", &data);
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
 }*/

/*std::vector<SpinWaveTwoBodyDecayTree> DecayGenerator::createAllValidTwoBodyDecayTrees() const {
 std::vector<SpinWaveTwoBodyDecayTree> final_decay_trees;

 // first create all seed decay trees
 std::vector<SpinWaveTwoBodyDecayTree> swtbdt = createSeedTwoBodyDecayTrees();

 std::vector<SpinWaveTwoBodyDecayTree> temp_swtbdt;

 while (swtbdt.size() > 0) {
 // then for all current seed decay trees, go through the complete list of
 // available decays and look for a possible decay
 // add it to the list and remove remaining particles entries
 for (unsigned int i = 0; i < swtbdt.size(); ++i) {
 if (swtbdt[i].available_particles_.size() > 1) {
 for (unsigned int available_particle_index = 1;
 available_particle_index < swtbdt[i].available_particles_.size();
 ++available_particle_index) {
 IndexList daughter_indices;
 daughter_indices.push_back(swtbdt[i].available_particles_[0]);
 daughter_indices.push_back(
 swtbdt[i].available_particles_[available_particle_index]);
 std::sort(daughter_indices.begin(), daughter_indices.end());
 auto result = two_body_decays_lookup_table_.find(daughter_indices);
 if (result != two_body_decays_lookup_table_.end()) {
 for (unsigned int twobodydecay_index = 0;
 twobodydecay_index < result->second.size();
 ++twobodydecay_index) {
 SpinWaveTwoBodyDecayTree temptree(swtbdt[i]);
 temptree.two_body_decays_.push_back(
 result->second[twobodydecay_index]);
 temptree.available_particles_.erase(
 temptree.available_particles_.begin()
 + available_particle_index);
 temptree.available_particles_.erase(
 temptree.available_particles_.begin());
 temptree.available_particles_.push_back(
 all_spin_wave_two_body_decays_[temptree.two_body_decays_.back()].mother_index_);


 temp_swtbdt.push_back(temptree);
 }
 }
 }
 }
 else {
 if (swtbdt[i].available_particles_.size() == 1) {
 final_decay_trees.push_back(swtbdt[i]);
 }
 }
 }
 swtbdt = temp_swtbdt;
 temp_swtbdt.clear();
 }      // do that until there are no remaining particles left

 final_decay_trees = filterDecayTreesForCorrectInitialState(final_decay_trees);

 return final_decay_trees;
 }

 std::vector<SpinWaveTwoBodyDecayTree> DecayGenerator::createSeedTwoBodyDecayTrees() const {
 std::vector<SpinWaveTwoBodyDecayTree> seed_decay_trees;
 std::vector<SpinWaveTwoBodyDecayTree> temp_seed_decay_trees;
 if (final_state_particles_.size() > 0) {

 for (unsigned int i = 0;
 i < final_state_particles_[0].spin_z_components_.size(); ++i) {
 ComPWA::ParticleProperties particle_properties =
 ComPWA::PhysConst::Instance().findParticle(
 final_state_particles_[0].particle_info_.name_);
 std::pair<bool, unsigned int> search_result = findSpinWaveIndex(
 particle_properties, final_state_particles_[0].spin_z_components_[i]);

 if (search_result.first) {
 unsigned int sw_index = search_result.second;
 SpinWaveTwoBodyDecayTree decay_tree;
 decay_tree.available_particles_.push_back(sw_index);
 seed_decay_trees.push_back(decay_tree);
 }
 }

 for (unsigned int fsp_index = 1; fsp_index < final_state_particles_.size();
 ++fsp_index) {
 for (unsigned int i = 0;
 i < final_state_particles_[fsp_index].spin_z_components_.size();
 ++i) {
 ComPWA::ParticleProperties particle_properties =
 ComPWA::PhysConst::Instance().findParticle(
 final_state_particles_[fsp_index].particle_info_.name_);

 std::pair<bool, unsigned int> search_result = findSpinWaveIndex(
 particle_properties,
 final_state_particles_[fsp_index].spin_z_components_[i]);

 if (search_result.first) {
 unsigned int sw_index = search_result.second;
 for (unsigned int existing_trees_index = 0;
 existing_trees_index < seed_decay_trees.size();
 ++existing_trees_index) {
 SpinWaveTwoBodyDecayTree decay_tree(
 seed_decay_trees[existing_trees_index]);
 decay_tree.available_particles_.push_back(sw_index);
 temp_seed_decay_trees.push_back(decay_tree);
 }
 }
 }
 seed_decay_trees = temp_seed_decay_trees;
 temp_seed_decay_trees.clear();
 }
 }
 return seed_decay_trees;
 }

 std::pair<bool, unsigned int> DecayGenerator::findSpinWaveIndex(
 const ComPWA::ParticleProperties& particle_properties,
 const Spin& spin_instance) const {
 SpinWave sw;
 sw.charge_ = particle_properties.charge_;
 sw.cparity_ = particle_properties.cparity_;
 sw.isospin_ = particle_properties.isospin_;
 sw.parity_ = particle_properties.parity_;
 sw.spin_ = spin_instance;

 // double check spin
 if (sw.spin_.J_numerator_ / sw.spin_.J_denominator_
 != particle_properties.spin_.J_numerator_
 / particle_properties.spin_.J_denominator_) {
 throw std::runtime_error(
 "DecayGenerator::findSpinWaveIndex: The specified final state wave spin does not equal to the one in the database!");
 }

 unsigned int position(0);
 bool found(false);

 auto result = std::find_if(all_spin_waves_.begin(), all_spin_waves_.end(),
 sw);

 if (result != all_spin_waves_.end()) {
 position = std::distance(all_spin_waves_.begin(), result);
 found = true;
 }
 else {
 std::stringstream ss;
 ss << "DecayGenerator::findSpinWaveIndex: The specified final state wave"
 << " cannot be found in\ the state pool please change your spin wave"
 << " pool generation constraints!" << std::endl;
 ss << sw;
 }

 return std::make_pair(found, position);
 }

 std::vector<SpinWaveTwoBodyDecayTree> DecayGenerator::filterDecayTreesForCorrectInitialState(
 const std::vector<SpinWaveTwoBodyDecayTree>& two_body_decay_trees) const {
 std::vector<SpinWaveTwoBodyDecayTree> filtered_decay_trees;

 IndexList mother_states;
 ComPWA::ParticleProperties particle_properties =
 ComPWA::PhysConst::Instance().findParticle(
 mother_state_particle_.particle_info_.name_);
 for (unsigned int i = 0; i < mother_state_particle_.spin_z_components_.size();
 ++i) {
 std::pair<bool, unsigned int> search_result = findSpinWaveIndex(
 particle_properties, mother_state_particle_.spin_z_components_[i]);
 if (search_result.first) {
 mother_states.push_back(search_result.second);
 }
 }

 for (unsigned int i = 0; i < two_body_decay_trees.size(); ++i) {
 bool good_decay_tree(false);
 if (two_body_decay_trees[i].available_particles_.size() == 1) {
 unsigned int spin_wave_index =
 two_body_decay_trees[i].available_particles_[0];
 for (unsigned int mother_state_index = 0;
 mother_state_index < mother_states.size(); ++mother_state_index) {
 if (spin_wave_index == mother_states[mother_state_index]) {
 good_decay_tree = true;
 break;
 }
 }
 }
 if (good_decay_tree)
 filtered_decay_trees.push_back(two_body_decay_trees[i]);
 }
 return filtered_decay_trees;
 }*/

void DecayGenerator::printDecayTree(const SpinWaveDecayTree& decay_tree) const {
  std::set<unsigned int> wave_list;
  for (auto decay = decay_tree.unique_decay_node_index_tree_.begin();
      decay != decay_tree.unique_decay_node_index_tree_.end(); ++decay) {
    std::cout << decay->first << " -> ";
    for (auto decay_daughter = decay->second.begin();
        decay_daughter != decay->second.end(); ++decay_daughter) {
      std::cout << *decay_daughter;
      if (decay_daughter != --decay->second.end()) {
        std::cout << " + ";
      }

    }
    std::cout << std::endl;
  }
  for (auto index_link =
      decay_tree.unique_decay_node_index_to_spin_wave_index_mapping_.begin();
      index_link
          != decay_tree.unique_decay_node_index_to_spin_wave_index_mapping_.end();
      ++index_link) {
    //std::cout << index_link->first << ": " << index_link->second
    //    << " but vector size is: " << all_spin_waves_.size() << std::endl;
    std::cout << index_link->first << ": "
        << all_spin_waves_[index_link->second] << std::endl;
  }
}

} /* namespace DecayTree */
} /* namespace Physics */
} /* namespace ComPWA */
