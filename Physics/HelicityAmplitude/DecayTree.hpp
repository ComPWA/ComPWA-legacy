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

#ifndef PHYSICS_HELICITYAMPLITUDE_DECAYTREE_HPP_
#define PHYSICS_HELICITYAMPLITUDE_DECAYTREE_HPP_

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/depth_first_search.hpp"
#include "boost/property_tree/ptree_fwd.hpp"

#include "ParticleStateDefinitions.hpp"

namespace HelicityFormalism {

struct DecayNode {
  ParticleStateInfo state_info_;
  boost::property_tree::ptree strength_and_phase_;

  bool operator==(const DecayNode &rhs) const {
    if (this->state_info_ != rhs.state_info_)
      return false;

    return true;
  }
  bool operator!=(const DecayNode &rhs) const {
    return !((*this) == rhs);
  }

  bool operator<(const DecayNode &rhs) const {
    if (this->state_info_ > rhs.state_info_)
      return false;
    else if (this->state_info_ < rhs.state_info_)
      return true;

    return true;
  }
  bool operator>(const DecayNode &rhs) const {
    return (rhs < *this);
  }
};

/**
 * Definition of our decay tree structure. We want a std::vector (vecS) as storage
 * for fast access times and low memory consumption. However then the uniqueness
 * of the nodes/vertices has to be ensured on our own (see addVertex() function).
 * The ParticleStates are used as the vertex property values.
 */
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
    DecayNode> HelicityTree;

class DecayTree {
  HelicityTree decay_tree_;

  std::vector<DecayNode> currently_grown_nodes_;

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

  struct ConnectedVertexDetector: public boost::dfs_visitor<> {
    ConnectedVertexDetector(unsigned int& connected) :
        connected_(connected) {
    }

    template<class Vertex, class Graph>
    void discover_vertex(Vertex v, Graph&) {
      ++connected_;
    }

  protected:
    unsigned int& connected_;
  };

  struct GetVertexDecendants: public boost::dfs_visitor<> {
    GetVertexDecendants(std::vector<unsigned int>& decendant_list,
        const boost::graph_traits<HelicityTree>::vertex_descriptor& vertex) :
        decendant_list_(decendant_list), vertex_(vertex), write_vertex_(false) {
    }

    template<class Vertex, class Graph>
    void discover_vertex(Vertex v, Graph& g) {
      //if (write_vertex_)
      //  decendant_list_.push_back(v);
      if (v == vertex_)
        write_vertex_ = true;
      if (write_vertex_)
        current_vertex_ = v;
    }

    template<class Vertex, class Graph>
    void finish_vertex(Vertex v, Graph& g) {
      if (write_vertex_ && v == current_vertex_)
        decendant_list_.push_back(v);
      if (v == vertex_)
        write_vertex_ = false;
    }

  protected:
    bool write_vertex_;
    boost::graph_traits<HelicityTree>::vertex_descriptor current_vertex_;
    const boost::graph_traits<HelicityTree>::vertex_descriptor& vertex_;
    std::vector<unsigned int>& decendant_list_;
  };

  struct GetAscendingVertexList: public boost::dfs_visitor<> {
    GetAscendingVertexList(
        std::vector<boost::graph_traits<HelicityTree>::vertex_descriptor>& vertices,
        bool with_leaves = true) :
        decendant_level(0), with_leaves_(with_leaves), vertices_(vertices) {
    }

    void finish_up() {
      std::list<boost::graph_traits<HelicityTree>::vertex_descriptor> sorted_list_of_vertices;
      for (auto iter = level_sorted_vertices_.begin();
          iter != level_sorted_vertices_.end(); ++iter) {
        sorted_list_of_vertices.insert(sorted_list_of_vertices.begin(),
            iter->second.begin(), iter->second.end());
      }
      vertices_.clear();
      vertices_.reserve(sorted_list_of_vertices.size());
      for (auto iter = sorted_list_of_vertices.begin();
          iter != sorted_list_of_vertices.end(); ++iter) {
        vertices_.push_back(*iter);
      }
    }

    template<class Vertex, class Graph>
    void discover_vertex(Vertex v, Graph& g) {
      ++decendant_level;
      increased_level = true;
    }

    template<class Vertex, class Graph>
    void finish_vertex(Vertex v, Graph& g) {
      --decendant_level;
      if (with_leaves_ || !increased_level) {
        level_sorted_vertices_[decendant_level].push_back(v);
      }
      increased_level = false;
      if (decendant_level == 0) {
        finish_up();
      }
    }

  protected:
    unsigned int decendant_level;
    bool increased_level;
    bool with_leaves_;
    std::vector<boost::graph_traits<HelicityTree>::vertex_descriptor>& vertices_;
    std::map<unsigned int,
        std::vector<boost::graph_traits<HelicityTree>::vertex_descriptor> > level_sorted_vertices_;
  };

  struct VertexWriter {
    VertexWriter(const HelicityTree& graph) :
        graph_(graph) {
    }
    void operator()(std::ostream& out,
        HelicityTree::vertex_descriptor v) const {
      out << "[label=\"" << graph_[v].state_info_.pid_information_.name_ << " ("
          << graph_[v].state_info_.spin_information_.J_ << ","
          << graph_[v].state_info_.spin_information_.M_ << ") \"]";
    }
  private:
    const HelicityTree& graph_;
  };

public:
  DecayTree();
  virtual ~DecayTree();

  bool isDisconnected() const;

  bool hasCycles() const;

  bool isDecayTreeValid() const;

  void clearCurrentGrownNodes();

  const HelicityTree& getHelicityDecayTree() const;

  void determineListOfDecayVertices();
  const std::vector<
      boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor>& getDecayVertexList() const;

  std::vector<
      boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor> getDecayNodesList() const;

  std::vector<DecayNode> getLowestLeaves() const;

  std::vector<DecayNode> getLeaves() const;

  boost::graph_traits<HelicityTree>::vertex_descriptor getTopNode() const;

  boost::graph_traits<HelicityTree>::vertex_descriptor addVertex(
      const DecayNode& node);

  void createDecay(const DecayNode &mother,
      const std::vector<ParticleStateInfo> &daughters);

  std::vector<ParticleStateInfo> createDecayProductsFinalStateParticleLists(
      const boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor& vertex) const;

  void print(std::ostream& os) const;

};

} /* namespace HelicityFormalism */

#endif /* PHYSICS_HELICITYAMPLITUDE_DECAYTREE_HPP_ */
