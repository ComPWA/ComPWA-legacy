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
namespace DecayTree {

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
  setupClipsEnvironment();

  createAllSpinWaves();
  createAllValidTwoBodyDecays();
  std::cout << "number of spin waves: " << all_spin_waves_.size() << std::endl;
  //std::cout << "number of two body decays: "
  //    << all_spin_wave_two_body_decays_.size() << std::endl;
  std::vector<SpinWaveTwoBodyDecayTree> two_body_decay_trees =
      createAllValidTwoBodyDecayTrees();

  std::cout << "we have " << two_body_decay_trees.size() << " decay trees!"
      << std::endl;

  for (unsigned int i = 0; i < two_body_decay_trees.size(); ++i) {
    printDecayTree(two_body_decay_trees[i]);
  }

  std::vector<ParticleIndexDecayTree> concrete_decay_trees;
  // now just port all the concrete realizations to an boost property tree
  // and write them to a xml file
  boost::property_tree::ptree decay_property_tree = portPropertyTree(
      concrete_decay_trees);
}

void DecayGenerator::createAllSpinWaves() {
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
  EnvRun(clips_environment, max_number_of_rules_to_fire);
  //EnvFacts(clips_environment, "stdout", NULL, -1, -1, -1);

  //high_resolution_clock::time_point t3 = high_resolution_clock::now();

  // get expertise
  std::vector<int> violated_rules_ids = getExpertise();

  //high_resolution_clock::time_point t4 = high_resolution_clock::now();

  EnvReset(clips_environment);

  //high_resolution_clock::time_point t5 = high_resolution_clock::now();

  /*duration<double> time_span1 = duration_cast<duration<double>>(t2 - t1);
   duration<double> time_span2 = duration_cast<duration<double>>(t3 - t2);
   duration<double> time_span3 = duration_cast<duration<double>>(t4 - t3);
   duration<double> time_span4 = duration_cast<duration<double>>(t5 - t4);

   std::cout<<"adding two body decay takes (seconds): "<< time_span1.count() <<std::endl;
   std::cout<<"running clips (seconds): "<< time_span2.count() <<std::endl;
   std::cout<<"getting results (seconds): "<< time_span3.count() <<std::endl;
   std::cout<<"resetting (seconds): "<< time_span4.count() <<std::endl;*/

  bool result(true);

  if (violated_rules_ids.size() > 0) {
    //std::cout << "this decay violates " << violated_rules_ids.size()
    //    << " rules\n";
    result = false;
  }

  return result;
}

void DecayGenerator::setupClipsEnvironment() const {
  std::stringstream path;
  path << getenv("COMPWA_DIR")
      << "/Physics/HelicityAmplitude/helicity_model.clp";
  EnvLoad(clips_environment, path.str().c_str());
  EnvReset(clips_environment);
}

void DecayGenerator::addTwoBodyDecayToClipsEnviroment(
    const SpinWaveTwoBodyDecay& two_body_decay) const {

  // add the SpinWaves
  addSpinWaveToClipsEnvironment(all_spin_waves_[two_body_decay.mother_index_]);
  addSpinWaveToClipsEnvironment(
      all_spin_waves_[two_body_decay.daughter_indices_.first]);
  addSpinWaveToClipsEnvironment(
      all_spin_waves_[two_body_decay.daughter_indices_.second]);

  void* two_body_decay_template;
  two_body_decay_template = EnvFindDeftemplate(clips_environment,
      "TwoBodyDecay");
  void* two_body_decay_fact = EnvCreateFact(clips_environment,
      two_body_decay_template);
  if (two_body_decay_fact != NULL) {
    DATA_OBJECT field;
    field.type = INTEGER;

    field.value = EnvAddLong(clips_environment,
        all_spin_waves_[two_body_decay.mother_index_].unique_id_);
    EnvPutFactSlot(clips_environment, two_body_decay_fact, "mother", &field);
    field.value = EnvAddLong(clips_environment,
        all_spin_waves_[two_body_decay.daughter_indices_.first].unique_id_);
    EnvPutFactSlot(clips_environment, two_body_decay_fact, "daughter1", &field);
    field.value = EnvAddLong(clips_environment,
        all_spin_waves_[two_body_decay.daughter_indices_.second].unique_id_);
    EnvPutFactSlot(clips_environment, two_body_decay_fact, "daughter2", &field);

    EnvAssert(clips_environment, two_body_decay_fact);
  }
}

void DecayGenerator::addSpinWaveToClipsEnvironment(
    const SpinWave& spin_wave) const {
  /*ComPWA::ParticleProperties particle_properties =
   ComPWA::PhysConst::Instance().findParticle(
   particle_info.particle_info_.name_);
   */
  void* spinwave_template;
  spinwave_template = EnvFindDeftemplate(clips_environment, "SpinWave");
// set the facts
  void* spinwave_fact = EnvCreateFact(clips_environment, spinwave_template);
  if (spinwave_fact != NULL) {
    DATA_OBJECT data_value;
    data_value.type = INTEGER;
    data_value.value = EnvAddLong(clips_environment, spin_wave.unique_id_);
    EnvPutFactSlot(clips_environment, spinwave_fact, "unique_id", &data_value);
    data_value.value = EnvAddLong(clips_environment, spin_wave.charge_);
    EnvPutFactSlot(clips_environment, spinwave_fact, "charge", &data_value);
    data_value.value = EnvAddLong(clips_environment,
        spin_wave.isospin_.J_numerator_);
    EnvPutFactSlot(clips_environment, spinwave_fact, "isospin_num",
        &data_value);
    data_value.value = EnvAddLong(clips_environment,
        spin_wave.isospin_.J_denominator_);
    EnvPutFactSlot(clips_environment, spinwave_fact, "isospin_denom",
        &data_value);
    data_value.value = EnvAddLong(clips_environment,
        spin_wave.isospin_.J_z_numerator_);
    EnvPutFactSlot(clips_environment, spinwave_fact, "isospin_z_num",
        &data_value);
    data_value.value = EnvAddLong(clips_environment,
        spin_wave.spin_.J_numerator_);
    EnvPutFactSlot(clips_environment, spinwave_fact, "spin_num", &data_value);
    data_value.value = EnvAddLong(clips_environment,
        spin_wave.spin_.J_denominator_);
    EnvPutFactSlot(clips_environment, spinwave_fact, "spin_denom", &data_value);
    data_value.value = EnvAddLong(clips_environment,
        spin_wave.spin_.J_z_numerator_);
    EnvPutFactSlot(clips_environment, spinwave_fact, "spin_z_num", &data_value);

    data_value.value = EnvAddLong(clips_environment, spin_wave.parity_);
    EnvPutFactSlot(clips_environment, spinwave_fact, "parity", &data_value);
    data_value.value = EnvAddLong(clips_environment, spin_wave.cparity_);
    EnvPutFactSlot(clips_environment, spinwave_fact, "cparity", &data_value);

    EnvAssignFactSlotDefaults(clips_environment, spinwave_fact);
    EnvAssert(clips_environment, spinwave_fact);
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
      clips_environment, clips_query.str().c_str(), &violated_rules_for_decay);

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

std::vector<SpinWaveTwoBodyDecayTree> DecayGenerator::createAllValidTwoBodyDecayTrees() const {
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

              /*decay_topology[]

              // if we only have 1 particle left this was the top node decay
              if(temptree.available_particles_.size() == 1) {

              }*/

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
}

void DecayGenerator::printDecayTree(
    const SpinWaveTwoBodyDecayTree& decay_tree) const {
  std::set<unsigned int> wave_list;
  for (unsigned int i = 0; i < decay_tree.two_body_decays_.size(); ++i) {
    const SpinWaveTwoBodyDecay& two_body_decay =
        all_spin_wave_two_body_decays_[decay_tree.two_body_decays_[i]];
    std::cout << two_body_decay.mother_index_ << " -> "
        << two_body_decay.daughter_indices_.first << " + "
        << two_body_decay.daughter_indices_.second << std::endl;
    wave_list.insert(two_body_decay.mother_index_);
    wave_list.insert(two_body_decay.daughter_indices_.first);
    wave_list.insert(two_body_decay.daughter_indices_.second);
  }
  for (auto iter = wave_list.begin(); iter != wave_list.end(); ++iter) {
    std::cout << *iter << std::endl;
    std::cout << all_spin_waves_[*iter] << std::endl;
  }
}

DecayConfiguration DecayGenerator::createDecayConfiguration(
    const std::vector<SpinWaveTwoBodyDecayTree>& two_body_decay_trees) const {
  DecayConfiguration decay_configuration;

  std::vector<ParticleIndexDecayTree> concrete_decay_trees;



  return decay_configuration;
}

boost::property_tree::ptree DecayGenerator::portPropertyTree(
    const std::vector<ParticleIndexDecayTree>& concrete_decay_trees) const {
  return boost::property_tree::ptree();
}

} /* namespace DecayTree */
} /* namespace ComPWA */
