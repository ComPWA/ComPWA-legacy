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

  const HelicityTree adjacency_list = decay_tree.getHelicityDecayTree();

  //typedef property_map<Graph, vertex_index_t>::type IndexMap;
  //    IndexMap index = get(vertex_index, g);

  std::pair<boost::graph_traits<HelicityTree>::vertex_iterator,
      boost::graph_traits<HelicityTree>::vertex_iterator> vp;

  std::pair<boost::graph_traits<HelicityTree>::out_edge_iterator,
      boost::graph_traits<HelicityTree>::out_edge_iterator> ep;

  /*for (vp = vertices(adjacency_list); vp.first != vp.second; ++vp.first) {

    ep = out_edges(vp.first, adjacency_list);

    HelicityAmplitude ha(adjacency_list[*vp.first],);

    amp_tree.amplitude_nodes.push_back(ha);
  }
*/
  return amp_tree;
}

} /* namespace HelicityFormalism */
