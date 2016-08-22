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
  //allowed_particle_names_.push_back("f0_980");
  //allowed_particle_names_.push_back("f0_1370");
  //allowed_particle_names_.push_back("f2_1270");
  allowed_particle_names_.push_back("omega");
  allowed_particle_names_.push_back("jpsi");
  //allowed_particle_names_.push_back("pi+");
  //allowed_particle_names_.push_back("pi-");
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
    auto particle_properties = ComPWA::PhysConst::Instance().findParticle(
        particle_info.particle_info_.name_);

    ComPWA::Spin s = particle_properties.getSpinLikeQuantumNumber(
        QuantumNumberIDs::SPIN);

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
    const SpinWaveDecayTree& two_body_decay_tree) {

  BOOST_LOG_TRIVIAL(debug)<<"adding spin wave decay tree to decay configuration \n";
  BOOST_LOG_TRIVIAL(debug)<<"number of decay nodes: "<< two_body_decay_tree.unique_decay_node_index_tree_.size() << std::endl;

  std::map<unsigned int, std::vector<ParticleStateInfo> > temp_state_pool;
  std::vector<
  std::vector<std::pair<ParticleStateInfo, std::vector<ParticleStateInfo> > > > decay_trees;

  for (auto const &decay_node : two_body_decay_tree.unique_decay_node_index_tree_) {
    BOOST_LOG_TRIVIAL(debug)<<"mother: "<<decay_node.first<<std::endl;
    if (temp_state_pool.find(decay_node.first) == temp_state_pool.end()) {
      temp_state_pool[decay_node.first] =
      createParticleStateInfoCandidates(
          two_body_decay_tree.unique_decay_node_index_to_spin_wave_index_mapping_.at(
              decay_node.first), decay_node.first);
      //if we could not find any particle for this mother state then just quit here
      if (temp_state_pool[decay_node.first].size() == 0)
      return;
    }

    std::vector<std::vector<ParticleStateInfo> > daughter_list_combinations;
    for (unsigned int i = 0; i < decay_node.second.size(); ++i) {
      BOOST_LOG_TRIVIAL(debug)<<"daugther: "<<decay_node.second[i]<<std::endl;
      if (temp_state_pool.find(decay_node.second[i]) == temp_state_pool.end()) {
        temp_state_pool[decay_node.second[i]] =
        createParticleStateInfoCandidates(
            two_body_decay_tree.unique_decay_node_index_to_spin_wave_index_mapping_.at(
                decay_node.second[i]), decay_node.second[i]);
        //if we could not find any particle for this daughter state then just quit here
        if (temp_state_pool[decay_node.second[i]].size() == 0)
        return;
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

    bool decay_trees_empty(false);
    if (decay_trees.size() == 0) {
      decay_trees_empty = true;
    }

    std::vector<
    std::vector<
    std::pair<ParticleStateInfo, std::vector<ParticleStateInfo> > > > new_decay_trees;

    BOOST_LOG_TRIVIAL(debug)<<"constructing decay trees and checking for single decay constraints...\n";
    // add them to the decay trees
    for (auto const& mother : temp_state_pool[decay_node.first]) {
      for (auto const& daughters : daughter_list_combinations) {
        // do some checks here
        if (!checkMass(mother, daughters))
        continue;

        BOOST_LOG_TRIVIAL(debug) << mother << std::endl;
        for (auto d : daughters) {
          BOOST_LOG_TRIVIAL(debug)<< d << std::endl;
        }

        if (decay_trees_empty) {
          std::vector<
          std::pair<ParticleStateInfo, std::vector<ParticleStateInfo> > > blank_decay_tree;
          auto decay_pair = std::make_pair(mother, daughters);
          blank_decay_tree.push_back(decay_pair);
          decay_trees.push_back(blank_decay_tree);
        }
        else {
          for (auto const& current_decay_tree : decay_trees) {
            auto decay_pair = std::make_pair(mother, daughters);
            if (isDecayValidForTree(decay_pair, current_decay_tree)) {
              std::vector<
              std::pair<ParticleStateInfo, std::vector<ParticleStateInfo> > > extended_decay_tree(
                  current_decay_tree);
              extended_decay_tree.push_back(decay_pair);
              new_decay_trees.push_back(extended_decay_tree);
            }
          }
        }
      }
    }
    if (!decay_trees_empty)
    decay_trees = new_decay_trees;
  }

  BOOST_LOG_TRIVIAL(debug)<<"checking for correct initial and final state and adding decay tree to decay configuration...\n";
  for (auto const& decay_tree : decay_trees) {
    if (!checkForCorrectIFState(decay_tree))
    continue;
    // reset particle picking pool
    current_total_particle_pool_ = total_particle_pool_;
    current_particle_mapping_.clear();
    std::cout << "decay tree\n";
    // for each possible mother state make copy for all pairs of daughters
    for (auto const &decay_node : decay_tree) {
      ParticleStateInfo mother = createParticleInstance(decay_node.first, isParticleIntermediateState(decay_node.first, decay_tree));
      BOOST_LOG_TRIVIAL(debug)<<" before!!!: trying to add decay for: "
      << decay_node.first.pid_information_.name_ << " ("
      << decay_node.first.unique_id_ << ") -> ";
      std::vector<ParticleStateInfo> daughters;
      for (auto const& daughter : decay_node.second) {
        daughters.push_back(createParticleInstance(daughter, isParticleIntermediateState(daughter, decay_tree)));
        BOOST_LOG_TRIVIAL(debug)<<daughter.pid_information_.name_ << " ("
        << daughter.unique_id_ << ") ";
      }
      BOOST_LOG_TRIVIAL(debug)<< std::endl;
      const boost::property_tree::ptree decay_strength_info_and_phase =
      createStrengthAndPhase();

      BOOST_LOG_TRIVIAL(debug)<< "trying to add decay for: " << mother.pid_information_.name_
      << " (" << mother.unique_id_ << ") -> ";
      for (auto daughter : daughters) {
        BOOST_LOG_TRIVIAL(debug)<< daughter.pid_information_.name_ << " ("
        << daughter.unique_id_ << ") ";
      }
      BOOST_LOG_TRIVIAL(debug)<< std::endl;

      decay_configuration.addDecayToCurrentDecayTree(mother, daughters,
          decay_strength_info_and_phase);
    }
    decay_configuration.addCurrentDecayTreeToList();
  }
  BOOST_LOG_TRIVIAL(debug)<<"finished converting spin wave decay trees to concrete decay trees...\n";
}

std::vector<ParticleStateInfo> DecayGenerator::createParticleStateInfoCandidates(
    unsigned int spin_wave_index, unsigned int unique_index) const {
  SpinWave spin_wave = all_spin_waves_[spin_wave_index];

  std::vector<ParticleStateInfo> candidates;

  // find candidates
  for (auto const& allowed_particle : allowed_particle_pool_) {
    if (spin_wave.supersetOf(allowed_particle)) {
      ParticleStateInfo ps;
      ps.spin_information_ = spin_wave.getSpinLikeQuantumNumber(
          QuantumNumberIDs::SPIN);

      ps.pid_information_.particle_id_ = allowed_particle.id_;
      ps.pid_information_.name_ = allowed_particle.name_;

      ps.dynamical_information_ = createDynamicInfo(allowed_particle,
          DynamicalFunctions::DynamicalInfoTypes::RELATIVE_BREIT_WIGNER);

      ps.unique_id_ = unique_index;

      candidates.push_back(ps);
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
    dynamical_info.put("type",
        ComPWA::Physics::DynamicalFunctions::DynamicalTypeToString.at(
            dynamical_type));
    dynamical_info.put("mass.value", particle_properties.mass_);
    dynamical_info.put("mass.fix", 1);
    dynamical_info.put("mass.min", 0.5 * particle_properties.mass_);
    dynamical_info.put("mass.max", 1.5 * particle_properties.mass_);
    dynamical_info.put("width.value", particle_properties.width_);
    dynamical_info.put("width.fix", 1);
    dynamical_info.put("width.min", 0.5 * particle_properties.width_);
    dynamical_info.put("width.max", 10.0 * particle_properties.width_);
    dynamical_info.put("mesonRadius", 1.0);
    dynamical_info.put("norm", 1);
    dynamical_info.put("par1", 1.0);
    dynamical_info.put("par2", 1.0);
  }

  return dynamical_info;
}

bool DecayGenerator::checkMass(const ParticleStateInfo& mother,
    const std::vector<ParticleStateInfo>& daughters) const {
  double sum_mass_daughters(0.0);
  for (auto daughter : daughters) {
    sum_mass_daughters += PhysConst::Instance().findParticle(
        daughter.pid_information_.particle_id_).mass_;
  }
  double mass_mother = PhysConst::Instance().findParticle(
      mother.pid_information_.particle_id_).mass_;
  return mass_mother > sum_mass_daughters;
}

bool DecayGenerator::isDecayValidForTree(
    const std::pair<ParticleStateInfo, std::vector<ParticleStateInfo> >& decay_pair,
    const std::vector<
        std::pair<ParticleStateInfo, std::vector<ParticleStateInfo> > >& decay_tree) const {
  // TODO: we need to check that the node we want to connect to is actually a free node, so a leaf
  // check if we can plug this decay on the bottom of our tree
  for (auto decay_node : decay_tree) {
    for (auto daughter : decay_node.second) {
      if (decay_pair.first == daughter)
        return true;
    }
  }
  // check if we can plug this decay on the top of our tree
  for (auto decay_node : decay_tree) {
    for (auto daughter : decay_pair.second) {
      if (decay_node.first == daughter)
        return true;
    }
  }
  return false;
}

bool DecayGenerator::checkForCorrectIFState(
    const std::vector<
        std::pair<ParticleStateInfo, std::vector<ParticleStateInfo> > >& decay_tree) const {
  ParticleStateInfo initial_state;
  std::vector<ParticleStateInfo> final_state;

  for (auto decay_node : decay_tree) {
    // check if this decay mother is initial state
    bool is_initial(true);
    for (auto other_decay_node : decay_tree) {
      for (auto other_daughter : other_decay_node.second) {
        if (other_daughter.unique_id_ == decay_node.first.unique_id_) {
          is_initial = false;
          break;
        }
      }
      if (!is_initial)
        break;
    }
    if (is_initial)
      initial_state = decay_node.first;

    // check if this decay daughers are final state ones
    for (auto daughter : decay_node.second) {
      bool is_final(true);
      for (auto other_decay_node : decay_tree) {
        if (daughter.unique_id_ == other_decay_node.first.unique_id_) {
          is_final = false;
          break;
        }
      }
      if (is_final)
        final_state.push_back(daughter);
    }
  }

  if (mother_state_particle_.particle_info_.particle_id_
      != initial_state.pid_information_.particle_id_)
    return false;

  std::vector<IFParticleInfo> true_final_state_copy(final_state_particles_);
  for (auto final_state_particle : final_state) {
    auto result =
        std::find_if(true_final_state_copy.begin(), true_final_state_copy.end(),
            [&](const IFParticleInfo& p) {return p.particle_info_.particle_id_ == final_state_particle.pid_information_.particle_id_;});
    if (result != true_final_state_copy.end()) {
      //remove this particle from pool and continue
      true_final_state_copy.erase(result);
    }
    else {
      return false;
    }
  }
  return true;
}

ParticleStateInfo DecayGenerator::createParticleInstance(
    const ParticleStateInfo& ps, bool make_coherent) {

  ParticleStateInfo return_ps(ps);

  auto used_result = current_particle_mapping_[ps.unique_id_].find(
      ps.pid_information_.particle_id_);
  if (used_result != current_particle_mapping_[ps.unique_id_].end()) {
    return_ps.unique_id_ = used_result->second.unique_id_;
  }
  else {
    auto result =
        std::find_if(current_total_particle_pool_.begin(),
            current_total_particle_pool_.end(),
            [&](const ParticleStateInfo& psi) {return psi.pid_information_.particle_id_ == ps.pid_information_.particle_id_;});

    // if we could not find anything then just add to the end off the total particle pool
    if (result == current_total_particle_pool_.end()) {
      total_particle_pool_.push_back(ps);
      total_particle_pool_[total_particle_pool_.size() - 1].unique_id_ =
          total_particle_pool_.size();
      return_ps.unique_id_ = total_particle_pool_[total_particle_pool_.size()
          - 1].unique_id_;
    }
    else {
      return_ps.unique_id_ = result->unique_id_;
      // remove this particle out of the current pool
      current_total_particle_pool_.erase(
          std::remove(current_total_particle_pool_.begin(),
              current_total_particle_pool_.end(), *result),
          current_total_particle_pool_.end());
    }
    current_particle_mapping_[ps.unique_id_][ps.pid_information_.particle_id_] =
        return_ps;
  }
  return_ps.coherent = make_coherent;

  return return_ps;
}

bool DecayGenerator::isParticleIntermediateState(const ParticleStateInfo& state,
    const std::vector<
        std::pair<ParticleStateInfo, std::vector<ParticleStateInfo> > >& decay_tree) const {
  // ok if its the top node or a leaf then it cant be intermediate
  bool is_leaf(true);
  bool is_top_node(true);
  for (auto const& decay_node : decay_tree) {
    if (decay_node.first == state) {
      is_leaf = false;
    }
    if (std::find(decay_node.second.begin(), decay_node.second.end(), state)
        != decay_node.second.end()) {
      is_top_node = false;
    }
  }
  return !(is_leaf || is_top_node);
}

const boost::property_tree::ptree DecayGenerator::createStrengthAndPhase() const {
  boost::property_tree::ptree decay_strength_info_and_phase;

  decay_strength_info_and_phase.put("strength.value", 1.0);
  decay_strength_info_and_phase.put("strength.fix", 0);
  decay_strength_info_and_phase.put("strength.min", 0);
  decay_strength_info_and_phase.put("strength.max", 5);
  decay_strength_info_and_phase.put("phase.value", 0.0);
  decay_strength_info_and_phase.put("phase.fix", 0);
  decay_strength_info_and_phase.put("phase.min", -100);
  decay_strength_info_and_phase.put("phase.max", 100);

  return decay_strength_info_and_phase;
}

std::vector<SpinWaveDecayTree> DecayGenerator::generate() {
  setupClipsEnvironment();

  addAllowedQuantumNumbersToClipsEnvironment();
  addConservedQuantumNumbersToClipsEnvironment();
  addAllInitialAndFinalStateCombinationsToClipsEnvironment();

// run
//EnvFacts(clips_environment_, "stdout", NULL, -1, -1, -1);
  int max_number_of_rules_to_fire(-1);
  EnvRun(clips_environment_, max_number_of_rules_to_fire);
// EnvFacts(clips_environment_, "stdout", NULL, -1, -1, -1);

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
  path << getenv("COMPWA_DIR") << "/Physics/DecayTree/Isospin.clp";
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
  dummy.required_quantum_numbers_names_ =
      allowed_qn.required_quantum_numbers_names_;
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
  SpinWave wave(
      ComPWA::PhysConst::Instance().findParticle(
          particle.particle_info_.name_));

  wave.spin_like_quantum_numbers_[QuantumNumberTranslator::Instance().getQuantumNumberName(
      ComPWA::QuantumNumberIDs::SPIN)] =
      particle.spin_z_components_[state_index];

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
    SpinWaveDecay sw_decay;

    DATA_OBJECT data;
    EnvGetFactSlot(clips_environment_, GetMFValue(GetValue(decay_fact_list), i),
        "mother", &data);
    unsigned int unique_id_mother = ValueToInteger(GetValue(data));
    sw_decay.mother_index_ =
        findSpinWaveListIndex(
            spin_wave_decay_tree.unique_decay_node_index_to_spin_wave_index_mapping_,
            unique_id_mother);

    IndexList daughter_ids;
    EnvGetFactSlot(clips_environment_, GetMFValue(GetValue(decay_fact_list), i),
        "daughters", &data);
    for (unsigned int daughter_index = GetDOBegin(data);
        daughter_index <= GetDOEnd(data); ++daughter_index) {
      unsigned int temp_daugher_index(
          ValueToInteger(GetMFValue(GetValue(data), daughter_index)));
      daughter_ids.push_back(temp_daugher_index);
      sw_decay.daughter_indices_.push_back(
          findSpinWaveListIndex(
              spin_wave_decay_tree.unique_decay_node_index_to_spin_wave_index_mapping_,
              temp_daugher_index));
    }

    EnvGetFactSlot(clips_environment_, GetMFValue(GetValue(decay_fact_list), i),
        "violating_quantum_number_list", &data);
    for (unsigned int violated_qn_index = GetDOBegin(data);
        violated_qn_index <= GetDOEnd(data); ++violated_qn_index) {
      sw_decay.violated_quantum_numbers_.push_back(
          QuantumNumberTranslator::Instance().getQuantumNumberEnum(
              ValueToString(GetMFValue(GetValue(data), violated_qn_index))));
    }
    spin_wave_decay_tree.unique_decay_node_index_tree_[unique_id_mother] =
        daughter_ids;

    unsigned int index_pos = std::find(all_spin_wave_decays_.begin(),
        all_spin_wave_decays_.end(), sw_decay) - all_spin_wave_decays_.begin();

    if (index_pos == all_spin_wave_decays_.size()) {
      index_pos = all_spin_wave_decays_.size();
      all_spin_wave_decays_.push_back(sw_decay);
    }

    spin_wave_decay_tree.unique_decay_indices_.push_back(index_pos);
    spin_wave_decay_tree.unique_decay_node_index_to_decay_index_mapping_[unique_id_mother] =
        spin_wave_decay_tree.unique_decay_indices_.size() - 1;
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
  auto result = std::find(all_spin_waves_.begin(), all_spin_waves_.end(), sw);
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
  spin.z_component_relevant = false;

  DATA_OBJECT spin_quantum_number_fact_list;
  std::stringstream clips_query;
  clips_query << "(find-fact ((?f SpinQuantumNumber)) (= ?f:unique_id "
      << spin_quantum_number_unique_id << "))";
  EnvEval(clips_environment_, clips_query.str().c_str(),
      &spin_quantum_number_fact_list);

  //EnvFacts(clips_environment_, "stdout", NULL, -1, -1, 50);

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

unsigned int DecayGenerator::findSpinWaveListIndex(
    const std::map<unsigned int, unsigned int>& index_mapping,
    unsigned int unique_index) const {
  auto result = index_mapping.find(unique_index);
  if (result != index_mapping.end())
    return result->second;
  else {
    std::stringstream ss;
    ss
        << "DecayGenerator::getDecayTreeFromClipsEnvironment: no spinwave found with index "
        << unique_index << "!";
    std::runtime_error(ss.str());
  }
}

void DecayGenerator::printDecayTree(const SpinWaveDecayTree& decay_tree) const {
  std::set<unsigned int> wave_list;
  for (auto const& decay : decay_tree.unique_decay_node_index_tree_) {
    unsigned int mother_index = decay.first;
    std::cout << mother_index << " -> ";
    for (auto decay_daughter = decay.second.begin();
        decay_daughter != decay.second.end(); ++decay_daughter) {
      std::cout << *decay_daughter;
      if (decay_daughter != --decay.second.end()) {
        std::cout << " + ";
      }

    }
    std::cout << std::endl;

    unsigned int decay_list_index(
        decay_tree.unique_decay_node_index_to_decay_index_mapping_.at(
            mother_index));
    auto sw_decay = all_spin_wave_decays_.at(
        decay_tree.unique_decay_indices_.at(decay_list_index));
    if (sw_decay.violated_quantum_numbers_.size() > 0) {
      std::cout << "violated quantum numbers: ";
      for (auto const& violated_qn : sw_decay.violated_quantum_numbers_) {
        std::cout
            << QuantumNumberTranslator::Instance().getQuantumNumberName(
                violated_qn);
      }
      std::cout << std::endl;
    }
  }
  for (auto const& index_link : decay_tree.unique_decay_node_index_to_spin_wave_index_mapping_) {
    //std::cout << index_link->first << ": " << index_link->second
    //    << " but vector size is: " << all_spin_waves_.size() << std::endl;
    std::cout << index_link.first << ": " << all_spin_waves_[index_link.second]
        << std::endl;
  }
}

} /* namespace DecayTree */
} /* namespace Physics */
} /* namespace ComPWA */
