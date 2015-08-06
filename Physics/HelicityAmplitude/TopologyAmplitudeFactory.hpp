#ifndef PHYSICS_HELICITYAMPLITUDE_TOPOLOGYAMPLITUDEFACTORY_HPP_
#define PHYSICS_HELICITYAMPLITUDE_TOPOLOGYAMPLITUDEFACTORY_HPP_

#include "TopologyAmplitude.hpp"
#include "DecayTree.hpp"

namespace HelicityFormalism {

class TopologyAmplitudeFactory {
  std::map<TwoBodyDecayInformation, std::shared_ptr<TwoBodyDecayAmplitude> > two_body_decay_amplitude_list_;
  std::map<unsigned int, std::shared_ptr<TwoBodyDecayAmplitude> > dynamical_function_list;    // TODO fix abs function stuff

  SequentialTwoBodyDecayAmplitude generateSequentialDecayAmplitude(
      const DecayTree& decay_tree);

  DecayTopology createDecayTopology(const DecayTree& decay_tree) const;

/*  std::vector<ParticleState> getConnectedFinalStateParticleListForVertex(
      const boost::graph_traits<HelicityTree>::vertex_descriptor& decay_vertex) const;

  void descendVertexAndFillConnectedFinalStateParticleList(
      const boost::graph_traits<HelicityTree>::vertex_descriptor& decay_vertex,
      std::vector<ParticleState>& final_state_particles) const;*/

public:
  TopologyAmplitudeFactory();
  virtual ~TopologyAmplitudeFactory();

  std::vector<TopologyAmplitude> generateTopologyAmplitudes(
      const std::vector<DecayTree>& decay_tree_collection);
};

} /* namespace HelicityFormalism */

#endif /* PHYSICS_HELICITYAMPLITUDE_TOPOLOGYAMPLITUDEFACTORY_HPP_ */
