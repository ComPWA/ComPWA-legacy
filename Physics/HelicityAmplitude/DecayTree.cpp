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

#include <boost/graph/connected_components.hpp>
#include <boost/graph/graphviz.hpp>

#include "DecayTree.hpp"

namespace HelicityFormalism {

DecayTree::DecayTree() {
}

DecayTree::~DecayTree() {
}

bool DecayTree::isDisconnected() const {
  unsigned int total_vertices = num_vertices(decay_tree_);
  std::vector<int> component(total_vertices);
  unsigned int connected_vertices = connected_components(decay_tree_,
      &component[0]);
  return total_vertices != connected_vertices;
}

bool DecayTree::hasCycles() const {
  bool has_cycle(false);
  CycleDetector vis(has_cycle);
  boost::depth_first_search(decay_tree_, boost::visitor(vis));
  return has_cycle;
}

bool DecayTree::isDecayTreeValid() const {
  if (hasCycles()) {
    throw std::runtime_error(
        "The decay tree has a cyclic dependency, meaning the tree is corrupted. Please fix the configuration file and rerun!");
  }
  if (isDisconnected()) {
    throw std::runtime_error(
        "The decay tree is not connected, meaning the tree is corrupted. Please fix the configuration file and rerun!");
  }
  return true;
}

void DecayTree::clearCurrentGrownNodes() {
  currently_grown_nodes_.clear();
}

const HelicityTree& DecayTree::getHelicityDecayTree() const {
  return decay_tree_;
}

std::vector<ParticleStateInfo> DecayTree::getLowestLeaves() const {
  return currently_grown_nodes_;
}

void DecayTree::determineListOfDecayVertices() {
  std::pair<boost::graph_traits<HelicityTree>::vertex_iterator,
      boost::graph_traits<HelicityTree>::vertex_iterator> vp;

  for (vp = vertices(decay_tree_); vp.first != vp.second; ++vp.first) {
    // if we have 1 or more edges going out of this vertex its a decay vertex
    if (0 < out_degree(*vp.first, decay_tree_))
      decay_vertex_list_.push_back(*vp.first);
  }
}

const std::vector<
    boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor>& DecayTree::getDecayVertexList() const {
  return decay_vertex_list_;
}

void DecayTree::createDecay(const ParticleStateInfo &mother,
    const std::vector<ParticleStateInfo> &daughters) {
// check if the particles already exist as vertices with addVertex()
// if so return vertex descriptor for this vertex, if not create a new
  boost::graph_traits<HelicityTree>::vertex_descriptor mother_vertex =
      addVertex(mother);
  // then make the correct inserts into the tree
  for (unsigned int i = 0; i < daughters.size(); ++i) {
    boost::add_edge(mother_vertex, addVertex(daughters[i]), decay_tree_);
    currently_grown_nodes_.push_back(daughters[i]);
  }
}

boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor DecayTree::addVertex(
    const ParticleStateInfo& particle) {
  boost::graph_traits<HelicityTree>::vertex_descriptor return_vertex;
  boost::graph_traits<HelicityTree>::vertex_iterator vi, vi_end, next;
  tie(vi, vi_end) = vertices(decay_tree_);
  for (next = vi; next != vi_end; ++next) {
    if (decay_tree_[*next] == particle) {
      return_vertex = *next;
      break;
    }
  }
  if (next == vi_end) {
    //particles_.push_back(particle);
    return_vertex = add_vertex(decay_tree_);
    decay_tree_[return_vertex] = particle;
  }
  return return_vertex;
}

std::vector<std::vector<ParticleStateInfo> > DecayTree::createDecayProductsFinalStateParticleLists(
    const boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor& vertex) const {

  std::vector<std::vector<ParticleStateInfo> > return_vector;

  std::pair<boost::graph_traits<HelicityTree>::out_edge_iterator,
      boost::graph_traits<HelicityTree>::out_edge_iterator> ep;

  ep = boost::out_edges(vertex, decay_tree_);

  unsigned final_state_particle_counter = 0;
  while (ep.first != ep.second) {
    std::vector<unsigned int> complete_decendant_list;

    GetVertexDecendants vis(complete_decendant_list, vertex);
    boost::depth_first_search(decay_tree_, boost::visitor(vis));

    std::vector<ParticleStateInfo> fs_particle_list;
    fs_particle_list.reserve(complete_decendant_list.size());
    for(unsigned int i = 0; i < complete_decendant_list.size(); ++i) {
      fs_particle_list.push_back(decay_tree_[complete_decendant_list[i]]);
    }

    std::sort(fs_particle_list.begin(), fs_particle_list.end());
    return_vector.push_back(fs_particle_list);
  }
  return return_vector;
}

void DecayTree::print(std::ostream& os) const {
// boost::associative_property_map<IndexNameMap> propmapIndex(index_label_map);
  VertexWriter vertex_writer(decay_tree_);
  write_graphviz(os, decay_tree_, vertex_writer);
}

} /* namespace HelicityFormalism */
