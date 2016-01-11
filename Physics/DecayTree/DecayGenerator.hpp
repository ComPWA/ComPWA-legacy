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

namespace ComPWA {

class ParticleProperties;

namespace DecayTree {

struct IFParticleInfo {
  unsigned int unique_id_;
  ComPWA::IDInfo particle_info_;
  std::vector<ComPWA::Spin> spin_z_components_;
};

struct SpinWave {
  unsigned int unique_id_;
  int charge_;
  Spin isospin_;
  Spin spin_;
  int parity_;
  int cparity_;

  bool operator()(const SpinWave& rhs) {
    if (charge_ != rhs.charge_)
      return false;
    if (isospin_ != rhs.isospin_)
      return false;
    if (spin_ != rhs.spin_)
      return false;
    if (parity_ != rhs.parity_)
      return false;
    if (cparity_ != rhs.cparity_)
      return false;
    return true;
  }

  friend std::ostream& operator<<(std::ostream& stream, const SpinWave& sw) {
    stream << "charge: " << sw.charge_ << std::endl;
    stream << "isospin: " << sw.isospin_.J_numerator_ << "/"
        << sw.isospin_.J_denominator_ << "(I_z= " << sw.isospin_.J_z_numerator_
        << "/" << sw.isospin_.J_denominator_ << ")" << std::endl;
    stream << "spin: " << sw.spin_.J_numerator_ << "/"
            << sw.spin_.J_denominator_ << "(J_z= " << sw.spin_.J_z_numerator_
            << "/" << sw.spin_.J_denominator_ << ")" << std::endl;
    stream << "c-parity: " << sw.cparity_ << std::endl;
    stream << "parity: " << sw.parity_ << std::endl;
    return stream;
  }
};

struct SpinWaveTwoBodyDecay {
  unsigned int mother_index_;
  IndexPair daughter_indices_;
};

struct SpinWaveTwoBodyDecayTree {
  unsigned int top_node_decay_index;
  std::map<unsigned int, IndexList> decay_index_mapping;
  std::vector<unsigned int> two_body_decays_;
  std::vector<unsigned int> available_particles_;
};

class DecayGenerator {
  friend class DecayGeneratorFacade;

  void *clips_environment;

  // containers and fields set by user
  std::vector<IFParticleInfo> final_state_particles_;
  IFParticleInfo mother_state_particle_;

  std::vector<ComPWA::Spin> allowed_spins_;
  std::vector<ComPWA::Spin> allowed_isospins_;
  std::vector<int> allowed_charges_;
  std::vector<int> allowed_parities_;
  std::vector<int> allowed_cparites_;

  // dynamically filled containers
  std::vector<ParticleStateInfo> total_particle_pool_;
  IndexList mother_state_particle_index_;
  IndexList final_state_particles_indices_;
  IndexList intermediate_state_particles_indices_;

  std::vector<SpinWave> all_spin_waves_;
  std::vector<SpinWaveTwoBodyDecay> all_spin_wave_two_body_decays_;
  std::map<IndexList, IndexList> two_body_decays_lookup_table_;

  void populateSpinZStates(IFParticleInfo& particle_info) const;

  int lookupPID(const std::string& name) const;

  void createAllSpinWaves();
  void createAllValidTwoBodyDecays();

  void setupClipsEnvironment() const;

  bool validateTwoBodyDecay(const SpinWaveTwoBodyDecay& two_body_decay);
  void addTwoBodyDecayToClipsEnviroment(
      const SpinWaveTwoBodyDecay& two_body_decay) const;
  void addSpinWaveToClipsEnvironment(const SpinWave& spin_wave) const;
  std::vector<int> getExpertise();
  ParticleStateInfo getSpinWaveFromClipsEnvironment(int unique_id) const;
  void setPIDInfo(ParticleStateInfo& particle, void* spin_wave_fact) const;
  void print(const ParticleStateInfo &psi) const;

  std::vector<SpinWaveTwoBodyDecayTree> createAllValidTwoBodyDecayTrees() const;
  std::vector<SpinWaveTwoBodyDecayTree> createSeedTwoBodyDecayTrees() const;
  std::pair<bool, unsigned int> findSpinWaveIndex(
      const ComPWA::ParticleProperties& particle_properties,
      const Spin& spin_instance) const;
  std::vector<SpinWaveTwoBodyDecayTree> filterDecayTreesForCorrectInitialState(
      const std::vector<SpinWaveTwoBodyDecayTree>& two_body_decay_trees) const;

  void printDecayTree(
      const SpinWaveTwoBodyDecayTree& decay_tree) const;

  DecayConfiguration createDecayConfiguration(
      const std::vector<SpinWaveTwoBodyDecayTree>& two_body_decay_trees) const;

  boost::property_tree::ptree portPropertyTree(
      const std::vector<ParticleIndexDecayTree>& concrete_decay_trees) const;

public:
  DecayGenerator();
  virtual ~DecayGenerator();

  IFParticleInfo createIFParticleInfo(const std::string& name) const;

  void addFinalStateParticles(const IFParticleInfo& particle_info);
  void setTopNodeState(const IFParticleInfo& particle_info);

  void generate();
};

} /* namespace DecayTree */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_DECAYGENERATOR_HPP_ */
