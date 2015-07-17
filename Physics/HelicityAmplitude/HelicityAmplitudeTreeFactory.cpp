#include "boost/graph/graph_traits.hpp"

#include "HelicityAmplitudeTreeFactory.hpp"
#include "HelicityDecayTree.hpp"

namespace HelicityFormalism {

HelicityAmplitudeTreeFactory::HelicityAmplitudeTreeFactory() {
  // TODO Auto-generated constructor stub

}

HelicityAmplitudeTreeFactory::~HelicityAmplitudeTreeFactory() {
  // TODO Auto-generated destructor stub
}

HelicityAmplitudeTree HelicityAmplitudeTreeFactory::generateAmplitudeTree(
    const HelicityDecayTree& decay_tree) const {
  HelicityAmplitudeTree amp_tree;

  const HelicityTree helicity_tree = decay_tree.getHelicityDecayTree();
  const std::vector<
      boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor>& decay_vertex_list =
      decay_tree.getDecayVertexList();

  const std::vector<
      boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor>::const_iterator decay_vertex_iter;

  for (decay_vertex_iter = decay_vertex_list.begin();
      decay_vertex_iter != decay_vertex_list.end(); ++decay_vertex_iter) {

    ParticleStatePair decay_products;

    std::pair<boost::graph_traits<HelicityTree>::out_edge_iterator,
        boost::graph_traits<HelicityTree>::out_edge_iterator> ep;

    ep = out_edges(*decay_vertex_iter, decay_tree_);

    unsigned final_state_particle_counter = 0;
    while (ep.first != ep.second) {
      ++final_state_particle_counter;
      boost::graph_traits<HelicityTree>::vertex_descriptor decay_product(
          target(ep.first, decay_tree_));
      if (1 == final_state_particle_counter) {
        decay_products.first = helicity_tree[decay_product];
      }
      else if (2 == final_state_particle_counter) {
        decay_products.second = helicity_tree[decay_product];
      }
      else {
        std::stringstream ss;
        ss
            << "HelicityAmplitudeTreeFactory: This decay vertex has more than 2 final state particles ("
            << final_state_particle_counter
            << " fs particles)! This is not supported in this model. Please fix the decay file.";
        throw std::runtime_error(ss.str());
      }
      ++ep.first;
    }

    HelicityAmplitude ha(helicity_tree[*decay_vertex_iter], decay_products);
    amp_tree.amplitude_nodes_.push_back(ha);
  }

  return amp_tree;
}

std::vector<std::vector<ParticleState> > HelicityAmplitudeTreeFactory::getConnectedFinalStateParticleListsForDecayVertices() const {
  std::vector<std::vector<ParticleState> > return_vector;

  const std::vector<
      boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor>::const_iterator decay_vertex_iter;

  for (decay_vertex_iter = decay_vertex_list.begin();
      decay_vertex_iter != decay_vertex_list.end(); ++decay_vertex_iter) {
    return_vector.push_back(
        getConnectedFinalStateParticleListForVertex(*decay_vertex_iter));
  }
  return return_vector;
}

std::vector<ParticleState> HelicityAmplitudeTreeFactory::getConnectedFinalStateParticleListForVertex(
    const boost::graph_traits<HelicityTree>::vertex_descriptor& decay_vertex) const {
  std::vector<ParticleState> final_state_particle_list;

  descendVertexAndFillConnectedFinalStateParticleList(decay_vertex,
      final_state_particle_list);
  return final_state_particle_list;
}

void HelicityAmplitudeTreeFactory::descendVertexAndFillConnectedFinalStateParticleList(
    const boost::graph_traits<HelicityTree>::vertex_descriptor& decay_vertex,
    std::vector<ParticleState>& final_state_particles) const {

  const HelicityTree helicity_tree = decay_tree.getHelicityDecayTree();

  std::pair<boost::graph_traits<HelicityTree>::out_edge_iterator,
      boost::graph_traits<HelicityTree>::out_edge_iterator> ep;

  for (ep = out_edges(decay_vertex, helicity_tree); ep.first != ep.second;
      ++ep.first) {
    boost::graph_traits<HelicityTree>::vertex_descriptor decay_product(
        target(ep.first, helicity_tree));
    if (0 == out_degree(decay_product, decay_tree_)) {
      final_state_particles.push_back(
          boost::graph_traits<HelicityTree>::vertex_descriptor);
    }
    else {
      descendVertexAndFillConnectedFinalStateParticleList(decay_product,
          final_state_particles);
    }
  }
}

} /* namespace HelicityFormalism */
