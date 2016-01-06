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

#include "Physics/HelicityAmplitude/ParticleStateDefinitions.hpp"
#include "Physics/HelicityAmplitude/DecayConfiguration.hpp"

namespace HelicityFormalism {

typedef std::pair<TwoBodyDecayTopology, std::vector<std::vector<IDInfo> > > TopologyRemainingStatePair;

struct IFParticleInfo {
  unsigned int unique_id_;
  IDInfo particle_info_;
  std::vector<ComPWA::Spin> spin_z_components_;
};

class DecayGenerator {
  void *clips_environment;

  // containers and fields set by user
  std::vector<IFParticleInfo> final_state_particles_;
  //std::vector<IDInfo> intermediate_state_particles_;
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

  void populateSpinZStates(IFParticleInfo& particle_info) const;

  //IndexList createParticleStates(const IFParticleInfo& particle_info);

  int lookupPID(const std::string& name) const;

  std::vector<TwoBodyDecayTopology> constructPossibleDecayTopologies() const;

  bool isTopologyBuildingFinished(
      const std::vector<TopologyRemainingStatePair>& current_topology_remaining_state_pairs) const;

  bool isTopologyNonExistant(
      const std::vector<TwoBodyDecayTopology> topology_pool,
      const TwoBodyDecayTopology& probe) const;

  std::vector<TopologyRemainingStatePair> constructNextLevel(
      const std::vector<TopologyRemainingStatePair>& current_topology_remaining_state_pairs) const;

  /* std::vector<ParticleIndexDecayTree> createDecayTreeRealizationTemplates(
   const std::vector<TwoBodyDecayTopology>& possible_decay_topologies) const;

   std::vector<TParticlePDG> createDecayNodeParticleCandiateList(
   const std::vector<IDInfo>& decay_products) const;

   std::vector<TParticlePDG> convertIDInfoToTParticlePDG(
   const std::vector<IDInfo>& particle_id_info_list) const;*/

  void createAllSpinWaves();
  void createAllTwoBodyDecays();
  void validateTwoBodyDecays();
  void createAllValidTwoBodyDecayTrees();

  void setupClipsEnvironment() const;

  void addIFParticleToClipsEnvironment(
      const IFParticleInfo& particle_info) const;

  void applyUserConditions() const;

  void getExpertise();

  ParticleStateInfo getSpinWaveFromClipsEnvironment(int unique_id) const;

  void setPIDInfo(ParticleStateInfo& particle, void* spin_wave_fact) const;

  void print(const ParticleStateInfo &psi) const;

  /*std::vector<ParticleIndexDecayTree> makeConcreteDecayTrees(
   const std::vector<ParticleIndexDecayTree>& decay_tree_realization_templates);*/

  boost::property_tree::ptree portPropertyTree(
      const std::vector<ParticleIndexDecayTree>& concrete_decay_trees) const;

public:
  DecayGenerator();
  virtual ~DecayGenerator();

  IFParticleInfo createIFParticleInfo(const std::string& name) const;

  void addFinalStateParticles(const IFParticleInfo& particle_info);
  //void addIntermediateStateParticles(const std::string& name);
  void setTopNodeState(const IFParticleInfo& particle_info);

  void generate();
};

} /* namespace HelicityFormalism */

#endif /* PHYSICS_HELICITYAMPLITUDE_DECAYGENERATOR_HPP_ */
