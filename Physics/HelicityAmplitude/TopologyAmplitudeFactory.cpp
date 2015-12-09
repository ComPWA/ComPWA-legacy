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
    name << helicity_tree[*decay_vertex_iter].state_info_.pid_information_.name_
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
            helicity_tree[decay_product].state_info_.pid_information_;
      }
      else if (2 == final_state_particle_counter) {
        decay_spin_info.final_state_.second =
            helicity_tree[decay_product].state_info_.spin_information_;
        decay_products_id_info.second =
            helicity_tree[decay_product].state_info_.pid_information_;
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
          two_body_decay_amplitude_list_[decay_spin_info] = std::shared_ptr
              < TwoBodyDecayAmplitude
              > (new TwoBodyDecayAmplitude(decay_spin_info));
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
                helicity_tree[*decay_vertex_iter].state_info_.pid_information_,
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
    return std::shared_ptr < DoubleParameter
        > (new DoubleParameter("fs_mass",
            PhysConst::Instance().findParticle(id_info.particle_id_).mass_));
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
      helicity_tree[decay_tree.getTopNode()].state_info_.pid_information_;

  std::vector<DecayNode> leaves = decay_tree.getLeaves();
  for (unsigned int i = 0; i < leaves.size(); ++i) {
    decay_topology.final_state_id_list_.push_back(
        leaves[i].state_info_.pid_information_);
  }

  for (decay_vertex_iter = decay_vertex_list.begin();
      decay_vertex_iter != decay_vertex_list.end(); ++decay_vertex_iter) {

    TwoBodyDecayIndices indices;

    boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor mother_vertex =
        findMotherVertex(*decay_vertex_iter, decay_tree);

    std::vector<ParticleStateInfo> fs_particle_list =
        decay_tree.createDecayProductsFinalStateParticleLists(mother_vertex);

    indices.mother_index_ = decay_topology.insertFinalStateContentList(
        convertToIDInfoList(fs_particle_list));

    fs_particle_list = decay_tree.createDecayProductsFinalStateParticleLists(
        *decay_vertex_iter);

    indices.decay_state_index_ = decay_topology.insertFinalStateContentList(
        convertToIDInfoList(fs_particle_list));

    decay_topology.final_state_content_unique_id_mapping_.insert(
        std::make_pair(helicity_tree[*decay_vertex_iter].state_info_.unique_id_,
            convertToUniqueIDList(fs_particle_list)));

    std::pair<boost::graph_traits<HelicityTree>::out_edge_iterator,
        boost::graph_traits<HelicityTree>::out_edge_iterator> ep;

    ep = boost::out_edges(*decay_vertex_iter, helicity_tree);

    IndexList daughters;

    unsigned final_state_particle_counter = 0;
    while (ep.first != ep.second) {
      fs_particle_list = decay_tree.createDecayProductsFinalStateParticleLists(
          boost::target(*ep.first, helicity_tree));

      daughters.push_back(
          helicity_tree[boost::target(*ep.first, helicity_tree)].state_info_.unique_id_);

      unsigned int index = decay_topology.insertFinalStateContentList(
          convertToIDInfoList(fs_particle_list));

      decay_topology.final_state_content_unique_id_mapping_.insert(
          std::make_pair(
              helicity_tree[boost::target(*ep.first, helicity_tree)].state_info_.unique_id_,
              convertToUniqueIDList(fs_particle_list)));

      if (final_state_particle_counter == 0) {
        indices.decay_products_.first = index;
      }
      else {
        indices.decay_products_.second = index;
      }
      ++final_state_particle_counter;
      ++ep.first;
    }

    decay_topology.decay_node_fs_content_index_infos_.push_back(indices);

    decay_topology.particle_unique_id_decay_tree_[helicity_tree[*decay_vertex_iter].state_info_.unique_id_] =
        daughters;

    if (indices.decay_state_index_ == indices.mother_index_)
      decay_topology.top_node_unique_id_ =
          helicity_tree[*decay_vertex_iter].state_info_.unique_id_;
  }

  // specify the evaluation order...
  // so the data has to be arranged in the same way the sequential amp is built
  // up
  const std::vector<
      boost::graph_traits<HelicityFormalism::HelicityTree>::vertex_descriptor>& sequential_decay_vertex_list =
      decay_tree.getDecayNodesList();

  for (decay_vertex_iter = sequential_decay_vertex_list.begin();
      decay_vertex_iter != sequential_decay_vertex_list.end();
      ++decay_vertex_iter) {
    std::pair<boost::graph_traits<HelicityTree>::out_edge_iterator,
        boost::graph_traits<HelicityTree>::out_edge_iterator> ep;

    ep = boost::out_edges(*decay_vertex_iter, helicity_tree);

    unsigned int final_state_particle_counter = 0;
    while (ep.first != ep.second) {
      ++final_state_particle_counter;
      ++ep.first;
    }
    if (final_state_particle_counter > 0) {
      if (helicity_tree[*decay_vertex_iter].state_info_.dynamical_information_.get<
          std::string>("type") != "topNode") {
        decay_topology.unique_id_decay_node_order_.push_back(
            helicity_tree[*decay_vertex_iter].state_info_.unique_id_);
      }
    }
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

std::vector<IDInfo> TopologyAmplitudeFactory::convertToIDInfoList(
    const std::vector<ParticleStateInfo>& fs_particle_list) const {
  std::vector<IDInfo> idinfo_list;
  for (auto iter = fs_particle_list.begin(); iter != fs_particle_list.end();
      ++iter) {
    idinfo_list.push_back(iter->pid_information_);
  }
  std::sort(idinfo_list.begin(), idinfo_list.end());
  return idinfo_list;
}

IndexList TopologyAmplitudeFactory::convertToUniqueIDList(
    const std::vector<ParticleStateInfo>& fs_particle_list) const {
  IndexList unique_id_list;
  for (auto iter = fs_particle_list.begin(); iter != fs_particle_list.end();
      ++iter) {
    unique_id_list.push_back(iter->unique_id_);
  }
  std::sort(unique_id_list.begin(), unique_id_list.end());
  return unique_id_list;
}

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
  p.pid = particle_state.pid_information_.particle_id_;
  return p;
}

} /* namespace HelicityFormalism */
