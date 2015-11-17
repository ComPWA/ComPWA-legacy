#include "boost/graph/graph_traits.hpp"

#include "Core/PhysConst.hpp"
#include "Physics/HelicityAmplitude/TopologyAmplitudeFactory.hpp"
#include "Physics/HelicityAmplitude/DecayTree.hpp"

namespace HelicityFormalism {

TopologyAmplitudeFactory::TopologyAmplitudeFactory() {
  // TODO Auto-generated constructor stub
}

TopologyAmplitudeFactory::~TopologyAmplitudeFactory() {
  // TODO Auto-generated destructor stub
}

HelicityFormalism::SequentialTwoBodyDecayAmplitude TopologyAmplitudeFactory::generateSequentialDecayAmplitude(
    const DecayTree& decay_tree) {
  SequentialTwoBodyDecayAmplitude full_decay_amplitude;

  const HelicityTree helicity_tree = decay_tree.getHelicityDecayTree();
  const std::vector<
      boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor>& decay_vertex_list =
      decay_tree.getDecayNodesList();

  std::vector<
      boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor>::const_iterator decay_vertex_iter;

  std::stringstream name;

// We need to build the amplitudes from the bottom up, because the resonances
// can have the mass of another daughter resonance as a parameter and we want
// to this parameter exists only once. Thats actually the reason the vertex
// list is presorted in ascending level order...
// Moreover we have a map storage for the resonance parameters, which have to
// be passed to the top level resonance as external parameters if needed
  std::map<IDInfo, ParameterList> resonance_parameter_lists;

  for (decay_vertex_iter = decay_vertex_list.begin();
      decay_vertex_iter != decay_vertex_list.end(); ++decay_vertex_iter) {
    // construct full decay tree amplitude name by intermediate state name concatenation
    name << helicity_tree[*decay_vertex_iter].state_info_.id_information_.name_
        << "_";

    TwoBodyDecaySpinInformation decay_spin_info;
    std::pair<IDInfo, IDInfo> decay_products_id_info;

    std::pair<boost::graph_traits<HelicityTree>::out_edge_iterator,
        boost::graph_traits<HelicityTree>::out_edge_iterator> ep;

    ep = boost::out_edges(*decay_vertex_iter, helicity_tree);

    unsigned int final_state_particle_counter = 0;
    while (ep.first != ep.second) {
      ++final_state_particle_counter;
      boost::graph_traits<HelicityTree>::vertex_descriptor decay_product(
          boost::target(*ep.first, helicity_tree));
      if (1 == final_state_particle_counter) {
        decay_spin_info.final_state_.first =
            helicity_tree[decay_product].state_info_.spin_information_;
        decay_products_id_info.first =
            helicity_tree[decay_product].state_info_.id_information_;
      }
      else if (2 == final_state_particle_counter) {
        decay_spin_info.final_state_.second =
            helicity_tree[decay_product].state_info_.spin_information_;
        decay_products_id_info.second =
            helicity_tree[decay_product].state_info_.id_information_;
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

    if (final_state_particle_counter > 0) {
      // TODO: top node is just ignored... and basically set as a constant factor 1.0
      if (helicity_tree[*decay_vertex_iter].state_info_.dynamical_information_.get<
          std::string>("type") != "topNode") {
        decay_spin_info.initial_state_ =
            helicity_tree[*decay_vertex_iter].state_info_.spin_information_;

        // create new amplitude if not yet existent
        if (two_body_decay_amplitude_list_.find(decay_spin_info)
            == two_body_decay_amplitude_list_.end()) {
          two_body_decay_amplitude_list_[decay_spin_info] = std::shared_ptr<
              TwoBodyDecayAmplitude>(
              new TwoBodyDecayAmplitude(decay_spin_info));
        }

        // do external parameter stuff here
        ParameterList daughter_masses_parameter_list;

        std::shared_ptr<DoubleParameter> daughter1_mass(
            new DoubleParameter("daughter1_mass",
                getResonanceMassParameter(resonance_parameter_lists,
                    decay_products_id_info.first)->GetValue()));
        std::shared_ptr<DoubleParameter> daughter2_mass(
            new DoubleParameter("daughter2_mass",
                getResonanceMassParameter(resonance_parameter_lists,
                    decay_products_id_info.second)->GetValue()));
        daughter_masses_parameter_list.AddParameter(daughter1_mass);
        daughter_masses_parameter_list.AddParameter(daughter2_mass);

        TwoBodyDecayInformation decay_info;
        decay_info.spin_info_ = decay_spin_info;
        decay_info.dynamical_info_.initial_state_ =
            helicity_tree[*decay_vertex_iter].state_info_.dynamical_information_;
        std::shared_ptr<DynamicalFunctions::AbstractDynamicalFunction> abs_function =
            dynamical_function_factory_.generateDynamicalFunction(decay_info,
                daughter_masses_parameter_list);
        full_decay_amplitude.full_decay_amplitude_chain_list_.push_back(
            std::make_pair(two_body_decay_amplitude_list_[decay_spin_info],
                abs_function));

        // create external parameter list for this decay, which can be
        // used by higher level decays..
        resonance_parameter_lists.insert(
            std::make_pair(
                helicity_tree[*decay_vertex_iter].state_info_.id_information_,
                abs_function->getParameterList()));

        // get strength and phase of this decay node
        const boost::property_tree::ptree& strength_and_phase_pt =
            helicity_tree[*decay_vertex_iter].strength_and_phase_;

        full_decay_amplitude.strength_ = generateDoubleParameter(
            strength_and_phase_pt.get_child("strength"), name.str());
        full_decay_amplitude.phase_ = generateDoubleParameter(
            strength_and_phase_pt.get_child("phase"), name.str());
      }
    }
  }

  return full_decay_amplitude;
}

std::shared_ptr<DoubleParameter> TopologyAmplitudeFactory::getResonanceMassParameter(
    const std::map<IDInfo, ParameterList>& resonance_parameter_lists,
    const IDInfo& id_info) const {
  // first check if we have an external parameter list for the daughters
  auto search_result = resonance_parameter_lists.find(id_info);
  if (search_result != resonance_parameter_lists.end()) {
    // if so, check for the resonance mass in that list
    return search_result->second.GetDoubleParameter("mass");
  }
  else {
    // if not then it has to be a final state particle so just create a
    // shared_ptr with that mass here
    PhysConst *physics_constants = PhysConst::instance();
    return std::shared_ptr<DoubleParameter>(
        new DoubleParameter("fs_mass",
            physics_constants->getMass(id_info.particle_id_)));
  }
}

std::shared_ptr<DoubleParameter> TopologyAmplitudeFactory::generateDoubleParameter(
    const boost::property_tree::ptree& pt, const std::string& name) const {
  std::shared_ptr<DoubleParameter> new_parameter(
      new DoubleParameter("mag_" + name, pt.get<double>("value"),
          pt.get<double>("min"), pt.get<double>("max")));
  new_parameter->FixParameter(pt.get<bool>("fix"));

  return new_parameter;
}

std::vector<TopologyAmplitude> TopologyAmplitudeFactory::generateTopologyAmplitudes(
    const std::vector<DecayTree>& decay_tree_collection) {
// first group decay trees according to their topology
  std::map<TwoBodyDecayTopology, std::vector<DecayTree> > topology_grouped_decay_trees;
  std::vector<DecayTree>::const_iterator decay_tree_iter;
  for (decay_tree_iter = decay_tree_collection.begin();
      decay_tree_iter != decay_tree_collection.end(); ++decay_tree_iter) {
    topology_grouped_decay_trees[createDecayTopology(*decay_tree_iter)].push_back(
        *decay_tree_iter);
  }

  std::vector<TopologyAmplitude> topology_amps;

  std::map<TwoBodyDecayTopology, std::vector<DecayTree> >::const_iterator grouped_decay_tree_iter;
  for (grouped_decay_tree_iter = topology_grouped_decay_trees.begin();
      grouped_decay_tree_iter != topology_grouped_decay_trees.end();
      ++grouped_decay_tree_iter) {
    TopologyAmplitude topology_amp;

    for (decay_tree_iter = grouped_decay_tree_iter->second.begin();
        decay_tree_iter != grouped_decay_tree_iter->second.end();
        ++decay_tree_iter) {
      topology_amp.sequential_decay_amplitude_list_.push_back(
          generateSequentialDecayAmplitude(*decay_tree_iter));
    }
    topology_amps.push_back(topology_amp);
  }

  return topology_amps;
}

std::vector<TwoBodyDecayTopology> TopologyAmplitudeFactory::generateDecayTopologies(
    std::vector<DecayTree>& decay_trees) const {
  std::vector<TwoBodyDecayTopology> decay_topologies;

  for (unsigned int i = 0; i < decay_trees.size(); ++i) {
    TwoBodyDecayTopology dt = createDecayTopology(decay_trees[i]);
    if (decay_topologies.end()
        == std::find(decay_topologies.begin(), decay_topologies.end(), dt)) {
      decay_topologies.push_back(dt);
    }
  }
  return decay_topologies;
}

TwoBodyDecayTopology TopologyAmplitudeFactory::createDecayTopology(
    const DecayTree& decay_tree) const {
  TwoBodyDecayTopology decay_topology;

  const HelicityTree helicity_tree = decay_tree.getHelicityDecayTree();
  const std::vector<
      boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor>& decay_vertex_list =
      decay_tree.getDecayVertexList();

  std::vector<
      boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor>::const_iterator decay_vertex_iter;

  decay_topology.top_node_id_info_ =
      helicity_tree[decay_tree.getTopNode()].state_info_.id_information_;

  std::vector<DecayNode> leaves = decay_tree.getLeaves();
  for (unsigned int i = 0; i < leaves.size(); ++i) {
    decay_topology.final_state_id_list_.push_back(
        leaves[i].state_info_.id_information_);
  }

  for (decay_vertex_iter = decay_vertex_list.begin();
      decay_vertex_iter != decay_vertex_list.end(); ++decay_vertex_iter) {

    TwoBodyDecayIndices indices;

    boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor mother_vertex =
        findMotherVertex(*decay_vertex_iter, decay_tree);

    std::vector<IDInfo> fs_particle_list =
        decay_tree.createDecayProductsFinalStateParticleLists(mother_vertex);

    indices.mother_index_ = decay_topology.insertFinalStateContentList(
        fs_particle_list);

    fs_particle_list = decay_tree.createDecayProductsFinalStateParticleLists(
        *decay_vertex_iter);

    indices.decay_state_index_ = decay_topology.insertFinalStateContentList(
        fs_particle_list);

    std::pair<boost::graph_traits<HelicityTree>::out_edge_iterator,
        boost::graph_traits<HelicityTree>::out_edge_iterator> ep;

    ep = boost::out_edges(*decay_vertex_iter, helicity_tree);

    unsigned final_state_particle_counter = 0;
    while (ep.first != ep.second) {
      fs_particle_list = decay_tree.createDecayProductsFinalStateParticleLists(
          boost::target(*ep.first, helicity_tree));

      unsigned int index = decay_topology.insertFinalStateContentList(
          fs_particle_list);

      if (final_state_particle_counter == 0) {
        indices.decay_products_.first = index;
      }
      else {
        indices.decay_products_.second = index;
      }
      ++final_state_particle_counter;
      ++ep.first;
    }

    decay_topology.decay_node_infos_.push_back(indices);
  }
  return decay_topology;
}

boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor TopologyAmplitudeFactory::findMotherVertex(
    const boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor & decay_node,
    const DecayTree& decay_tree) const {
  // try to find mother of this decay_vertex
  boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor mother(
      decay_node);

  const HelicityTree helicity_tree = decay_tree.getHelicityDecayTree();
  const std::vector<
      boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor>& decay_vertex_list =
      decay_tree.getDecayVertexList();
  // loop over vertices

  std::vector<
      boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor>::const_iterator decay_vertex_iter;

  for (decay_vertex_iter = decay_vertex_list.begin();
      decay_vertex_iter != decay_vertex_list.end(); ++decay_vertex_iter) {

    std::pair<boost::graph_traits<HelicityTree>::out_edge_iterator,
        boost::graph_traits<HelicityTree>::out_edge_iterator> ep;

    ep = boost::out_edges(*decay_vertex_iter, helicity_tree);

    unsigned final_state_particle_counter = 0;
    while (ep.first != ep.second) {
      if (helicity_tree[decay_node].state_info_
          == helicity_tree[boost::target(*ep.first, helicity_tree)].state_info_) {
        return *decay_vertex_iter;
      }
      ++ep.first;
    }
  }
  return mother;
}

/*TwoBodyDecayTopology TopologyAmplitudeFactory::createDecayTopology(
 const DecayTree& decay_tree) const {
 TwoBodyDecayTopology decay_topology;

 const HelicityTree helicity_tree = decay_tree.getHelicityDecayTree();
 const std::vector<
 boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor>& decay_vertex_list =
 decay_tree.getDecayVertexList();

 std::vector<
 boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor>::const_iterator decay_vertex_iter;

 std::cout << "number of decay vertices for decay tree "
 << decay_vertex_list.size() << std::endl;
 decay_tree.print(std::cout);

 decay_topology.top_node_id_info_ =
 helicity_tree[decay_tree.getTopNode()].state_info_.id_information_;

 std::vector<DecayNode> leaves = decay_tree.getLeaves();
 for (unsigned int i = 0; i < leaves.size(); ++i) {
 decay_topology.final_state_id_list_.push_back(
 leaves[i].state_info_.id_information_);
 }

 for (decay_vertex_iter = decay_vertex_list.begin();
 decay_vertex_iter != decay_vertex_list.end(); ++decay_vertex_iter) {

 std::vector<std::vector<IDInfo> > fs_particle_lists =
 decay_tree.createDecayProductsFinalStateParticleLists(
 *decay_vertex_iter);

 TwoBodyDecayIndices indices;
 std::vector<IDInfo> mother_fs_particle_list;

 std::sort(fs_particle_lists.begin(), fs_particle_lists.end());
 // if (fs_particle_lists.size() == 2) {
 for (unsigned int i = 0; i < fs_particle_lists.size(); ++i) {
 std::vector<std::vector<IDInfo> >::const_iterator result = std::find(
 decay_topology.final_state_content_id_lists_.begin(),
 decay_topology.final_state_content_id_lists_.end(),
 fs_particle_lists[i]);

 unsigned int index = result
 - decay_topology.final_state_content_id_lists_.begin();
 if (result == decay_topology.final_state_content_id_lists_.end()) {
 decay_topology.final_state_content_id_lists_.push_back(
 fs_particle_lists[i]);
 }

 if (i == 0) {
 indices.decay_products_.first = index;
 }
 else {
 indices.decay_products_.second = index;
 }

 mother_fs_particle_list.insert(mother_fs_particle_list.end(),
 fs_particle_lists[i].begin(), fs_particle_lists[i].end());
 }
 /* }
 else if (fs_particle_lists.size() != 0) {
 std::stringstream ss;
 ss
 << "HelicityAmplitudeTreeFactory: This decay vertex does not have 2 final state particles ("
 << fs_particle_lists.size()
 << " fs particles)! This is not supported in this model. Please fix the decay file.";
 throw std::runtime_error(ss.str());
 }*/

// do the same for the mother
/* std::vector<std::vector<IDInfo> >::const_iterator result = std::find(
 decay_topology.final_state_content_id_lists_.begin(),
 decay_topology.final_state_content_id_lists_.end(),
 mother_fs_particle_list);

 unsigned int mother_index = result
 - decay_topology.final_state_content_id_lists_.begin();
 if (result == decay_topology.final_state_content_id_lists_.end()) {
 decay_topology.final_state_content_id_lists_.push_back(
 mother_fs_particle_list);
 }
 indices.mother_index_ = mother_index;

 decay_topology.decay_node_infos_.push_back(indices);

 decay_topology.print();
 }
 return decay_topology;
 }*/

Event TopologyAmplitudeFactory::createDummyEvent(
    const DecayTree& decay_tree) const {
  Event dummy_event;
  std::vector<DecayNode> final_state_particles = decay_tree.getLeaves();
  for (unsigned int i = 0; i < final_state_particles.size(); ++i) {
    dummy_event.addParticle(
        createParticle(final_state_particles[i].state_info_));
  }
  return dummy_event;
}

Particle TopologyAmplitudeFactory::createParticle(
    const ParticleStateInfo& particle_state) const {
  Particle p;
  p.pid = particle_state.id_information_.particle_id_;
  return p;
}

} /* namespace HelicityFormalism */
