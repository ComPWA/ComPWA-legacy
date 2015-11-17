//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

#ifndef PHYSICS_HELICITYAMPLITUDE_TOPOLOGYAMPLITUDEFACTORY_HPP_
#define PHYSICS_HELICITYAMPLITUDE_TOPOLOGYAMPLITUDEFACTORY_HPP_

#include "Core/Particle.hpp"
#include "Core/Event.hpp"

#include "Physics/HelicityAmplitude/TopologyAmplitude.hpp"
#include "Physics/HelicityAmplitude/DecayTree.hpp"
#include "Physics/DynamicalDecayFunctions/DynamicalFunctionFactory.hpp"

namespace HelicityFormalism {

class TopologyAmplitudeFactory {
  DynamicalFunctions::DynamicalFunctionFactory dynamical_function_factory_;

  std::map<TwoBodyDecaySpinInformation, std::shared_ptr<TwoBodyDecayAmplitude> > two_body_decay_amplitude_list_;

  HelicityFormalism::SequentialTwoBodyDecayAmplitude generateSequentialDecayAmplitude(
      const DecayTree& decay_tree);

  std::shared_ptr<DoubleParameter> getResonanceMassParameter(
      const std::map<IDInfo, ParameterList>& resonance_parameter_lists,
      const IDInfo& id_info) const;

  std::shared_ptr<DoubleParameter> generateDoubleParameter(
      const boost::property_tree::ptree& pt, const std::string& name) const;

  TwoBodyDecayTopology createDecayTopology(const DecayTree& decay_tree) const;

  std::pair<std::vector<ParticleStateInfo>, std::vector<ParticleStateInfo> > createDecayProductsFinalStateParticleLists(
      const boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor& vertex) const;

  boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor findMotherVertex(
      const boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor & decay_node,
      const DecayTree& graph) const;

  /*  std::vector<ParticleState> getConnectedFinalStateParticleListForVertex(
   const boost::graph_traits<HelicityTree>::vertex_descriptor& decay_vertex) const;

   void descendVertexAndFillConnectedFinalStateParticleList(
   const boost::graph_traits<HelicityTree>::vertex_descriptor& decay_vertex,
   std::vector<ParticleState>& final_state_particles) const;*/

  Particle createParticle(const ParticleStateInfo& particle_state) const;

public:
  TopologyAmplitudeFactory();
  virtual ~TopologyAmplitudeFactory();

  std::vector<TopologyAmplitude> generateTopologyAmplitudes(
      const std::vector<DecayTree>& decay_tree_collection);

  std::vector<TwoBodyDecayTopology> generateDecayTopologies(
      std::vector<DecayTree>& decay_trees) const;

  Event createDummyEvent(const DecayTree& decay_tree) const;
};

} /* namespace HelicityFormalism */

#endif /* PHYSICS_HELICITYAMPLITUDE_TOPOLOGYAMPLITUDEFACTORY_HPP_ */
