//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//		Stefan Pflueger - initial implementation
//-------------------------------------------------------------------------------

#ifndef PHYSICS_HELICITYAMPLITUDE_DECAYGENERATOR_HPP_
#define PHYSICS_HELICITYAMPLITUDE_DECAYGENERATOR_HPP_

#include <vector>

#include "Core/Utility.hpp"

#include "Physics/DecayTree/DecayConfiguration.hpp"

#include "Physics/DynamicalDecayFunctions/DynamicalFunctionFactory.hpp"

namespace ComPWA {

class ParticleProperties;

namespace Physics {
namespace DecayTree {

enum class QuantumNumberTypes {
  SINGLE_PARTICLE_BASED, COMPOSITE_PARTICLE_BASED
};

template<class T> struct AllowedQuantumNumbers {
  std::string quantum_number_name_;
  std::vector<T> allowed_values_;
  std::vector<std::string> required_quantum_numbers_names_;
  QuantumNumberTypes type;
};

struct SpinWave {
  std::map<std::string, ComPWA::Spin> spin_like_quantum_numbers_;
  std::map<std::string, int> integer_like_quantum_numbers_;
  std::map<std::string, double> double_like_quantum_numbers_;

  bool operator()(const SpinWave& rhs) {
    if (spin_like_quantum_numbers_ != rhs.spin_like_quantum_numbers_)
      return false;
    if (integer_like_quantum_numbers_ != rhs.integer_like_quantum_numbers_)
      return false;
    if (double_like_quantum_numbers_ != rhs.double_like_quantum_numbers_)
      return false;
    return true;
  }

  friend std::ostream& operator<<(std::ostream& stream, const SpinWave& sw) {
    for (auto iter = sw.spin_like_quantum_numbers_.begin();
        iter != sw.spin_like_quantum_numbers_.end(); ++iter) {
      stream << iter->first << ": " << iter->second.J_numerator_ << "/"
          << iter->second.J_denominator_ << " (z="
          << iter->second.J_z_numerator_ << ")\n";
    }
    for (auto iter = sw.integer_like_quantum_numbers_.begin();
        iter != sw.integer_like_quantum_numbers_.end(); ++iter) {
      stream << iter->first << ": " << iter->second << "\n";
    }
    return stream;
  }
};

struct IFParticleInfo {
  unsigned int unique_id_;
  ComPWA::IDInfo particle_info_;
  std::vector<ComPWA::Spin> spin_z_components_;
};

struct SpinWaveDecay {
  unsigned int mother_index_;
  IndexList daughter_indices_;
  std::vector<ComPWA::QuantumNumbers> violated_quantum_numbers_;

  bool operator==(const SpinWaveDecay& rhs) const {
    if (mother_index_
        != rhs.mother_index_)
      return false;
    if (daughter_indices_ != rhs.daughter_indices_)
      return false;
    if (violated_quantum_numbers_
        != rhs.violated_quantum_numbers_)
      return false;
    return true;
  }
};

struct SpinWaveDecayTree {
  unsigned int top_node_unique_decay_node_index_;
  std::map<unsigned, IndexList> unique_decay_node_index_tree_;
  std::map<unsigned int, unsigned int> unique_decay_node_index_to_spin_wave_index_mapping_;
  std::vector<unsigned int> unique_decay_indices_;
  std::map<unsigned int, unsigned int> unique_decay_node_index_to_decay_index_mapping_;

  bool operator==(const SpinWaveDecayTree& rhs) const {
    if (top_node_unique_decay_node_index_
        != rhs.top_node_unique_decay_node_index_)
      return false;
    if (unique_decay_node_index_tree_ != rhs.unique_decay_node_index_tree_)
      return false;
    if (unique_decay_node_index_to_spin_wave_index_mapping_
        != rhs.unique_decay_node_index_to_spin_wave_index_mapping_)
      return false;
    return true;
  }
};

class DecayGenerator {
  friend class DecayGeneratorFacade;

  void *clips_environment_;

  // containers and fields set by user
  std::vector<IFParticleInfo> final_state_particles_;
  IFParticleInfo mother_state_particle_;

  unsigned int unique_spin_state_counter_;
  std::vector<AllowedQuantumNumbers<ComPWA::Spin> > allowed_spin_like_quantum_numbers_;
  std::vector<AllowedQuantumNumbers<int> > allowed_integer_like_quantum_numbers_;

  std::vector<std::string> conserved_quantum_numbers_;

  std::vector<std::string> allowed_particle_names_;

  // dynamically filled containers
  std::vector<ParticleStateInfo> total_particle_pool_;
  std::vector<ParticleProperties> allowed_particle_pool_;
  //IndexList mother_state_particle_index_;
  //IndexList final_state_particles_indices_;
  //IndexList intermediate_state_particles_indices_;

  std::vector<SpinWave> all_spin_waves_;
  std::vector<SpinWaveDecay> all_spin_wave_decays_;
  std::map<IndexList, IndexList> two_body_decays_lookup_table_;

  void populateSpinZStates(IFParticleInfo& particle_info) const;

  int lookupPID(const std::string& name) const;

  //void createAllSpinWaves();
  //void createAllValidTwoBodyDecays();

  void setupClipsEnvironment() const;

  void addAllowedQuantumNumbersToClipsEnvironment();
  void addSpinLikeAllowedQuantumNumbersToClipsEnvironment(
      const AllowedQuantumNumbers<ComPWA::Spin>& allowed_qn);
  unsigned int addSpinQuantumNumberToClipsEnvironment(
      const ComPWA::Spin& spin_state);
  void addIntLikeAllowedQuantumNumbersToClipsEnvironment(
      const AllowedQuantumNumbers<int>& allowed_qn) const;

  void addConservedQuantumNumbersToClipsEnvironment();

  void addAllInitialAndFinalStateCombinationsToClipsEnvironment();
  SpinWave createSpinWave(const IFParticleInfo& particle,
      unsigned int state_index) const;
  void addInitialAndFinalStateToClipsEnvironment(const SpinWave& initial_state,
      const std::vector<SpinWave>& final_state);
  void* addSpinWaveToClipsEnvironment(const SpinWave& spinwave);

  std::vector<SpinWaveDecayTree> getExpertise();
  SpinWaveDecayTree getDecayTreeFromClipsEnvironment(void* decay_tree_fact);
  unsigned int createSpinWave(void* spin_wave_fact);
  ComPWA::Spin createSpin(int spin_quantum_number_unique_id) const;
  unsigned int findSpinWaveListIndex(
      const std::map<unsigned int, unsigned int>& index_mapping,
      unsigned int unique_index) const;

  void initializeAllowedParticlePool();

  void addSpinWaveTwoBodyDecayToDecayConfiguration(
      DecayConfiguration& decay_configuration,
      const SpinWaveDecayTree& two_body_decay_tree) const;

  std::vector<ParticleStateInfo> createParticleStateInfoCandidates(
      unsigned int spin_wave_index) const;

  DynamicalInfo createDynamicInfo(const ParticleProperties& particle_properties,
      ComPWA::Physics::DynamicalFunctions::DynamicalInfoTypes dynamical_type) const;

  /* bool validateTwoBodyDecay(const SpinWaveTwoBodyDecay& two_body_decay);
   void addTwoBodyDecayToClipsEnviroment(
   const SpinWaveTwoBodyDecay& two_body_decay) const;
   void addSpinWaveToClipsEnvironment(const SpinWave& spin_wave) const;
   std::vector<int> getExpertise();
   ParticleStateInfo getSpinWaveFromClipsEnvironment(int unique_id) const;
   void setPIDInfo(ParticleStateInfo& particle, void* spin_wave_fact) const;
   void print(const ParticleStateInfo &psi) const;*/

  /* std::vector<SpinWaveTwoBodyDecayTree> createAllValidTwoBodyDecayTrees() const;
   std::vector<SpinWaveTwoBodyDecayTree> createSeedTwoBodyDecayTrees() const;
   std::pair<bool, unsigned int> findSpinWaveIndex(
   const ComPWA::ParticleProperties& particle_properties,
   const Spin& spin_instance) const;
   std::vector<SpinWaveTwoBodyDecayTree> filterDecayTreesForCorrectInitialState(
   const std::vector<SpinWaveTwoBodyDecayTree>& two_body_decay_trees) const;
   */
  void printDecayTree(const SpinWaveDecayTree& decay_tree) const;

  std::vector<SpinWaveDecayTree> generate();

public:
  DecayGenerator();
  virtual ~DecayGenerator();

  IFParticleInfo createIFParticleInfo(const std::string& name) const;

  void addFinalStateParticles(const IFParticleInfo& particle_info);
  void setTopNodeState(const IFParticleInfo& particle_info);

  DecayConfiguration createDecayConfiguration();
};

} /* namespace DecayTree */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_DECAYGENERATOR_HPP_ */
