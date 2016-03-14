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

#include "boost/graph/connected_components.hpp"
#include "boost/graph/graphviz.hpp"
#include "boost/property_tree/ptree.hpp"

#include "DecayTree.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayTree {

DecayTree::DecayTree() {
}

DecayTree::~DecayTree() {
}

bool DecayTree::isDisconnected() const {
  unsigned int total_vertices = num_vertices(decay_tree_);
  unsigned int connected_vertices(0);
  ConnectedVertexDetector vis(connected_vertices);
  boost::depth_first_search(decay_tree_, boost::visitor(vis));
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

std::vector<DecayNode> DecayTree::getLowestLeaves() const {
  return currently_grown_nodes_;
}

std::vector<DecayNode> DecayTree::getLeaves() const {
  std::vector<unsigned int> complete_decendant_list;
  GetVertexDecendants vis(complete_decendant_list, getTopNode());
  boost::depth_first_search(decay_tree_, boost::visitor(vis));

  std::vector<DecayNode> fs_particle_list;
  fs_particle_list.reserve(complete_decendant_list.size());
  for (unsigned int i = 0; i < complete_decendant_list.size(); ++i) {
    fs_particle_list.push_back(decay_tree_[complete_decendant_list[i]]);
  }

  return fs_particle_list;
}

boost::graph_traits<HelicityTree>::vertex_descriptor DecayTree::getTopNode() const {
  boost::graph_traits<HelicityTree>::vertex_descriptor top_node;
  for (auto vertex = decay_vertex_list_.begin();
      vertex != decay_vertex_list_.end(); ++vertex) {
    // go through edges and check that this vertex never occurred in any edge as target
    std::pair<boost::graph_traits<HelicityTree>::edge_iterator,
        boost::graph_traits<HelicityTree>::edge_iterator> ep = boost::edges(
        decay_tree_);

    bool not_top_node(false);

    for (auto edge = ep.first; edge != ep.second; ++edge) {
      if (boost::target(*edge, decay_tree_) == *vertex) {
        not_top_node = true;
        break;
      }
    }
    if (!not_top_node) {
      top_node = *vertex;
      break;
    }
  }
  return top_node;
}

void DecayTree::determineListOfDecayVertices() {
  decay_vertex_list_.clear();

  GetAscendingVertexList vis(decay_vertex_list_, false);
  boost::depth_first_search(decay_tree_, boost::visitor(vis));
}

const std::vector<
    boost::graph_traits<HelicityTree>::vertex_descriptor>& DecayTree::getDecayVertexList() const {
  return decay_vertex_list_;
}

std::vector<
    boost::graph_traits<HelicityTree>::vertex_descriptor> DecayTree::getDecayNodesList() const {
  std::vector<boost::graph_traits<HelicityTree>::vertex_descriptor> decay_nodes_list;
  GetAscendingVertexList vis(decay_nodes_list);
  boost::depth_first_search(decay_tree_, boost::visitor(vis));
  return decay_nodes_list;
}

void DecayTree::createDecay(const DecayNode &mother,
    const std::vector<ParticleStateInfo> &daughters) {
// check if the particles already exist as vertices with addVertex()
// if so return vertex descriptor for this vertex, if not create a new
  boost::graph_traits<HelicityTree>::vertex_descriptor mother_vertex =
      addVertex(mother);
// then make the correct inserts into the tree
  for (unsigned int i = 0; i < daughters.size(); ++i) {
    DecayNode daughter;
    daughter.state_info_ = daughters[i];
    boost::add_edge(mother_vertex, addVertex(daughter), decay_tree_);
    currently_grown_nodes_.push_back(daughter);
  }

  determineListOfDecayVertices();
}

boost::graph_traits<HelicityTree>::vertex_descriptor DecayTree::addVertex(
    const DecayNode& node) {
  boost::graph_traits<HelicityTree>::vertex_descriptor return_vertex;
  boost::graph_traits<HelicityTree>::vertex_iterator vi, vi_end, next;
  tie(vi, vi_end) = vertices(decay_tree_);
  for (next = vi; next != vi_end; ++next) {
    if (decay_tree_[*next] == node) {
      // if we found thats node already
      // then update strength and phase info if available
      if (!node.strength_and_phase_.empty()) {
        decay_tree_[*next].strength_and_phase_ = node.strength_and_phase_;
      }
      return_vertex = *next;
      break;
    }
  }
  if (next == vi_end) {
    //particles_.push_back(particle);
    return_vertex = add_vertex(decay_tree_);
    decay_tree_[return_vertex] = node;
  }
  return return_vertex;
}

std::vector<ParticleStateInfo> DecayTree::createDecayProductsFinalStateParticleLists(
    const boost::graph_traits<HelicityTree>::vertex_descriptor& vertex) const {
  std::vector<unsigned int> complete_decendant_list;

  GetVertexDecendants vis(complete_decendant_list, vertex);
  boost::depth_first_search(decay_tree_, boost::visitor(vis));

  std::vector<ParticleStateInfo> fs_particle_list;
  fs_particle_list.reserve(complete_decendant_list.size());
  for (unsigned int i = 0; i < complete_decendant_list.size(); ++i) {
    fs_particle_list.push_back(
        decay_tree_[complete_decendant_list[i]].state_info_);
  }

  return fs_particle_list;
}

void DecayTree::print(std::ostream& os) const {
// boost::associative_property_map<IndexNameMap> propmapIndex(index_label_map);
  VertexWriter vertex_writer(decay_tree_);
  write_graphviz(os, decay_tree_, vertex_writer);
}

} /* namespace DecayTree */
} /* namespace Physics */
} /* namespace ComPWA */
