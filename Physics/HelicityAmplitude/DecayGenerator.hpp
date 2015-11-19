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

#include "Physics/HelicityAmplitude/ParticleStateDefinitions.hpp"
#include "Physics/HelicityAmplitude/DecayConfiguration.hpp"

namespace HelicityFormalism {

typedef std::pair<TwoBodyDecayTopology, std::vector<std::vector<IDInfo> > > TopologyRemainingStatePair;

class DecayGenerator {
  // containers and fields set by user
  std::vector<IDInfo> final_state_particles_;
  std::vector<IDInfo> intermediate_state_particles_;
  IDInfo mother_state_particle_;
  std::vector<Spin> allowed_L_;
  std::vector<Spin> allowed_M_;

  // dynamically filled containers
  std::vector<IDInfo> total_particle_pool_;



  IDInfo createIDInfo(const std::string& name);

  int lookupPID(const std::string& name) const;

  std::vector<TwoBodyDecayTopology> constructPossibleDecayTopologies() const;

  bool isTopologyBuildingFinished(
      const std::vector<TopologyRemainingStatePair>& current_topology_remaining_state_pairs) const;

  bool isTopologyNonExistant(
      const std::vector<TwoBodyDecayTopology> topology_pool,
      const TwoBodyDecayTopology& probe) const;

  std::vector<TopologyRemainingStatePair> constructNextLevel(
      const std::vector<TopologyRemainingStatePair>& current_topology_remaining_state_pairs) const;

  std::vector<ParticleIndexDecayTree> createDecayTreeRealizationTemplates(
      const std::vector<TwoBodyDecayTopology>& possible_decay_topologies) const;

  std::vector<ParticleIndexDecayTree> makeConcreteDecayTrees(
      const std::vector<ParticleIndexDecayTree>& decay_tree_realization_templates);

  boost::property_tree::ptree portPropertyTree(
      const std::vector<ParticleIndexDecayTree>& concrete_decay_trees) const;

public:
  DecayGenerator();
  virtual ~DecayGenerator();

  void addFinalStateParticles(const std::string& name);
  void addIntermediateStateParticles(const std::string& name);
  void setTopNodeState(const std::string& name);

  void generate();
};

} /* namespace HelicityFormalism */

#endif /* PHYSICS_HELICITYAMPLITUDE_DECAYGENERATOR_HPP_ */
