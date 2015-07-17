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

#ifndef HELICITYDECAYTREE_HPP_
#define HELICITYDECAYTREE_HPP_

#include "HelicityStateDefinitions.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>

namespace HelicityFormalism {

/**
 * Definition of our decay tree structure. We want a std::vector (vecS) as storage
 * for fast access times and low memory consumption. However then the uniqueness
 * of the nodes/vertices has to be ensured on our own (see addVertex() function).
 * The ParticleStates are used as the vertex property values.
 */
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
    ParticleState> HelicityTree;

class HelicityDecayTree {
  HelicityTree decay_tree_;

  std::vector<ParticleState> currently_grown_nodes_;

  boost::graph_traits<HelicityTree>::vertex_descriptor current_vertex_;

  boost::graph_traits<HelicityTree>::vertex_descriptor top_level_vertex_;
  std::vector<boost::graph_traits<HelicityTree>::vertex_descriptor> decay_vertex_list_;

  struct CycleDetector: public boost::dfs_visitor<> {
    CycleDetector(bool& has_cycle) :
        has_cycle_(has_cycle) {
    }

    template<class Vertex, class Graph>
    void discover_vertex(Vertex v, Graph&) {
      vertex_counter_[v] = 1;
    }

    template<class Edge, class Graph>
    void examine_edge(Edge e, Graph& g) {
      vertex_counter_[target(e, g)]++;
      if (vertex_counter_[target(e, g)] > 1)
        has_cycle_ = true;
    }

    template<class Edge, class Graph>
    void back_edge(Edge, Graph&) {
      has_cycle_ = true;
    }

  protected:
    std::map<unsigned int, unsigned int> vertex_counter_;
    bool& has_cycle_;
  };

  /*struct FinalStateParticleFinder: public boost::dfs_visitor<> {
   FinalStateParticleFinder() {
   }

   template<class Vertex, class Graph>
   void finish_vertex(v, g) {
   final_state_vertices_.push_back(v);
   }

   std::vector<boost::graph_traits<HelicityTree>::vertex_descriptor> final_state_vertices_;
   };*/

  struct VertexWriter {
    VertexWriter(const HelicityTree& graph) :
        graph_(graph) {
    }
    void operator()(std::ostream& out,
        HelicityTree::vertex_descriptor v) const {
      out << "[label=\"" << graph_[v].name_ << "\"]";
    }
  private:
    const HelicityTree& graph_;
  };

public:
  HelicityDecayTree();
  virtual ~HelicityDecayTree();

  bool isDisconnected() const;

  bool hasCycles() const;

  bool isDecayTreeValid() const;

  void clearCurrentGrownNodes();

  const HelicityTree& getHelicityDecayTree() const;

  void determineListOfDecayVertices();
  std::vector<
      boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor>& getDecayVertexList() const;

  std::vector<ParticleState> getLowestLeaves() const;

  boost::graph_traits<HelicityTree>::vertex_descriptor addVertex(
      const ParticleState& particle);

  void createDecay(const ParticleState &mother,
      const ParticleStatePair &daughters);

  void print(std::ostream& os) const;
};

} /* namespace HelicityFormalism */

#endif /* HELICITYDECAYTREE_HPP_ */
