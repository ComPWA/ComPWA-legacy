#include "boost/graph/graph_traits.hpp"

#include "Core/PhysConst.hpp"

#include "Physics/HelicityAmplitude/TopologyAmplitudeFactory.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

using ComPWA::Physics::DecayTree::DecayTree;
using ComPWA::Physics::DecayTree::DecayNode;
using ComPWA::Physics::DecayTree::HelicityTree;

TopologyAmplitudeFactory::TopologyAmplitudeFactory() :
    two_body_decay_amplitude_list_(), dynamical_function_factory_() {
  // TODO Auto-generated constructor stub
}

TopologyAmplitudeFactory::~TopologyAmplitudeFactory() {
  // TODO Auto-generated destructor stub
}

std::vector<TopologyAmplitude> TopologyAmplitudeFactory::generateTopologyAmplitudes(
    const std::vector<DecayTree>& decay_tree_collection) {
  // clear old parameter list
  global_parameter_list_.clear();

  // first group decay trees according to their topology
  std::vector<TwoBodyDecayTopology> decay_topologies = generateDecayTopologies(
      decay_tree_collection);
  std::vector<std::vector<DecayTree> > topology_grouped_decay_trees(
      decay_topologies.size());

  for (auto const& decay_tree : decay_tree_collection) {
    auto result = std::find(decay_topologies.begin(), decay_topologies.end(),
        createDecayTopology(decay_tree));

    unsigned int position = result - decay_topologies.begin();

    topology_grouped_decay_trees[position].push_back(decay_tree);
  }

  std::vector<TopologyAmplitude> topology_amps;

  for (auto const& grouped_decay_tree : topology_grouped_decay_trees) {
    TopologyAmplitude topology_amp;
    for (auto const& decay_tree : grouped_decay_tree) {
      topology_amp.sequential_decay_amplitude_list_.push_back(
          generateSequentialDecayAmplitude(decay_tree));
    }
    topology_amps.push_back(topology_amp);
  }

  return topology_amps;
}

std::vector<TwoBodyDecayTopology> TopologyAmplitudeFactory::generateDecayTopologies(
    const std::vector<DecayTree>& decay_trees) const {
  std::vector<TwoBodyDecayTopology> decay_topologies;

  for (auto const& decay_tree : decay_trees) {
    TwoBodyDecayTopology dt = createDecayTopology(decay_tree);
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
  const std::vector<boost::graph_traits<HelicityTree>::vertex_descriptor>& decay_vertex_list =
      decay_tree.getDecayVertexList();

  std::vector<boost::graph_traits<HelicityTree>::vertex_descriptor>::const_iterator decay_vertex_iter;

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

    boost::graph_traits<HelicityTree>::vertex_descriptor mother_vertex =
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
  const std::vector<boost::graph_traits<HelicityTree>::vertex_descriptor>& sequential_decay_vertex_list =
      decay_tree.getDecayNodesList();

  for (decay_vertex_iter = sequential_decay_vertex_list.begin();
      decay_vertex_iter != sequential_decay_vertex_list.end();
      ++decay_vertex_iter) {
    /* std::pair<boost::graph_traits<HelicityTree>::out_edge_iterator,
     boost::graph_traits<HelicityTree>::out_edge_iterator> ep;

     ep = boost::out_edges(*decay_vertex_iter, helicity_tree);

     unsigned int final_state_particle_counter = 0;
     while (ep.first != ep.second) {
     ++final_state_particle_counter;
     ++ep.first;
     }*/

    if (!decay_tree.isDecayVertexALeaf(*decay_vertex_iter)) {
      //if (helicity_tree[*decay_vertex_iter].state_info_.dynamical_information_.get<
      //    std::string>("type") != "topNode") {
      decay_topology.unique_id_decay_node_order_.push_back(
          helicity_tree[*decay_vertex_iter].state_info_.unique_id_);
      //}
    }
  }

  return decay_topology;
}

boost::graph_traits<HelicityTree>::vertex_descriptor TopologyAmplitudeFactory::findMotherVertex(
    const boost::graph_traits<HelicityTree>::vertex_descriptor & decay_node,
    const DecayTree& decay_tree) const {
// try to find mother of this decay_vertex
  boost::graph_traits<HelicityTree>::vertex_descriptor mother(decay_node);

  const HelicityTree helicity_tree = decay_tree.getHelicityDecayTree();
  const std::vector<boost::graph_traits<HelicityTree>::vertex_descriptor>& decay_vertex_list =
      decay_tree.getDecayVertexList();
// loop over vertices

  std::vector<boost::graph_traits<HelicityTree>::vertex_descriptor>::const_iterator decay_vertex_iter;

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

SequentialTwoBodyDecayAmplitude TopologyAmplitudeFactory::generateSequentialDecayAmplitude(
    const DecayTree& decay_tree) {
  SequentialTwoBodyDecayAmplitude full_decay_amplitude;

  const HelicityTree helicity_tree = decay_tree.getHelicityDecayTree();
  const std::vector<boost::graph_traits<HelicityTree>::vertex_descriptor>& decay_vertex_list =
      decay_tree.getDecayNodesList();

  std::vector<boost::graph_traits<HelicityTree>::vertex_descriptor>::const_iterator decay_vertex_iter;

// We need to build the amplitudes from the bottom up, because the resonances
// can have the mass of another daughter resonance as a parameter and we want
// to this parameter exists only once. Thats actually the reason the vertex
// list is presorted in ascending level order...

  for (decay_vertex_iter = decay_vertex_list.begin();
      decay_vertex_iter != decay_vertex_list.end(); ++decay_vertex_iter) {
    // construct full decay tree amplitude name by intermediate state name concatenation
    std::stringstream name;

    appendParticleInfoToName(name,
        helicity_tree[*decay_vertex_iter].state_info_);

    std::stringstream parity_related_name(name.str());

    TwoBodyDecaySpinInformation decay_spin_info;
    std::pair<ParticleStateInfo, ParticleStateInfo> decay_products_ps_info;

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
        decay_products_ps_info.first = helicity_tree[decay_product].state_info_;
      }
      else if (2 == final_state_particle_counter) {
        decay_spin_info.final_state_.second =
            helicity_tree[decay_product].state_info_.spin_information_;
        decay_products_ps_info.second =
            helicity_tree[decay_product].state_info_;
      }
      else {
        std::stringstream ss;
        ss
            << "HelicityAmplitudeTreeFactory: This decay vertex has more than 2 final state particles ("
            << final_state_particle_counter
            << " fs particles)! This is not supported in this model. Please fix the decay file.";
        throw std::runtime_error(ss.str());
      }
      appendParticleInfoToName(name, helicity_tree[decay_product].state_info_);

      ParticleStateInfo temp_psi(helicity_tree[decay_product].state_info_);
      temp_psi.spin_information_.J_z_numerator_ =
          -temp_psi.spin_information_.J_z_numerator_;
      appendParticleInfoToName(parity_related_name,
          helicity_tree[decay_product].state_info_);

      ++ep.first;
    }

    if (final_state_particle_counter > 0) {
      decay_spin_info.initial_state_ =
          helicity_tree[*decay_vertex_iter].state_info_.spin_information_;

      FullTwoBodyDecayAmplitude full_two_body_decay_amplitude;

      // create new amplitude if not yet existent
      if (two_body_decay_amplitude_list_.find(decay_spin_info)
          == two_body_decay_amplitude_list_.end()) {
        two_body_decay_amplitude_list_[decay_spin_info] = std::shared_ptr<
            TwoBodyDecayAmplitude>(new TwoBodyDecayAmplitude(decay_spin_info));
      }
      full_two_body_decay_amplitude.angular_part_ =
          two_body_decay_amplitude_list_[decay_spin_info];

      // do external parameter stuff here
      DynamicalFunctions::ExternalParameters daughter_masses_parameter_list;

      if (helicity_tree[*decay_vertex_iter].state_info_.dynamical_information_.get<
          std::string>("type") != "topNode") {
        fixParticlePropertyTree(
            decay_products_ps_info.first.dynamical_information_,
            decay_products_ps_info.first.pid_information_);
        std::shared_ptr<DoubleParameter> daughter1_mass(
            generateGlobalParameter(
                decay_products_ps_info.first.dynamical_information_.get_child(
                    "mass"),
                decay_products_ps_info.first.pid_information_.name_ + "_mass",
                "").first);

        fixParticlePropertyTree(
            decay_products_ps_info.second.dynamical_information_,
            decay_products_ps_info.second.pid_information_);
        std::shared_ptr<DoubleParameter> daughter2_mass(
            generateGlobalParameter(
                decay_products_ps_info.second.dynamical_information_.get_child(
                    "mass"),
                decay_products_ps_info.second.pid_information_.name_ + "_mass",
                "").first);

        daughter_masses_parameter_list.parameters_.AddParameter(daughter1_mass);
        daughter_masses_parameter_list.parameter_name_mapping_["daughter1_mass"] =
            daughter1_mass->GetName();

        daughter_masses_parameter_list.parameters_.AddParameter(daughter2_mass);
        daughter_masses_parameter_list.parameter_name_mapping_["daughter2_mass"] =
            daughter2_mass->GetName();

      }

      std::shared_ptr<DynamicalFunctions::AbstractDynamicalFunction> abs_function =
          dynamical_function_factory_.generateDynamicalFunction(
              helicity_tree[*decay_vertex_iter].state_info_,
              daughter_masses_parameter_list);

      full_two_body_decay_amplitude.dynamical_part_ = abs_function;

      // create external parameter list for this decay, which can be
      // used by higher level decays..
      auto param_list = abs_function->getParameterList();
      for (auto const& param : param_list.GetDoubleParameters())
        global_parameter_list_.insert(std::make_pair(param->GetName(), param));

      // get strength and phase of this decay node
      const boost::property_tree::ptree &strength_and_phase_pt =
          helicity_tree[*decay_vertex_iter].strength_and_phase_;

      auto strength_param = generateGlobalParameter(
          strength_and_phase_pt.get_child("strength"), name.str() + "_mag",
          parity_related_name.str() + "_mag");

      auto phase_param = generateGlobalParameter(
          strength_and_phase_pt.get_child("phase"), name.str() + "_phase",
          parity_related_name.str() + "_phase");

      full_two_body_decay_amplitude.strength_ = strength_param.first;
      full_two_body_decay_amplitude.phase_ = phase_param.first;

      if (strength_param.second) {
        full_decay_amplitude.factor *=
            getParityFactor(
                helicity_tree[*decay_vertex_iter].state_info_.pid_information_.particle_id_,
                decay_products_ps_info.first.pid_information_.particle_id_,
                decay_products_ps_info.second.pid_information_.particle_id_);
      }

      full_decay_amplitude.decay_amplitude.push_back(
          full_two_body_decay_amplitude);
    }
  }

  return full_decay_amplitude;
}

double TopologyAmplitudeFactory::getParityFactor(int pid_mother, int pid_d1,
    int pid_d2) const {
  int parity_product(1);

  auto const& particle1 = PhysConst::Instance().findParticle(pid_d1);
  parity_product *= particle1.getIntLikeQuantumNumber(QuantumNumberIDs::PARITY);
  auto spin = particle1.getSpinLikeQuantumNumber(QuantumNumberIDs::SPIN);
  double spin_sum(spin.J_numerator_ / spin.J_denominator_);

  auto const& particle2 = PhysConst::Instance().findParticle(pid_d2);
  parity_product *= particle2.getIntLikeQuantumNumber(QuantumNumberIDs::PARITY);
  spin = particle2.getSpinLikeQuantumNumber(QuantumNumberIDs::SPIN);
  spin_sum += spin.J_numerator_ / spin.J_denominator_;

  auto const& particle = PhysConst::Instance().findParticle(pid_mother);
  parity_product *= particle.getIntLikeQuantumNumber(QuantumNumberIDs::PARITY);
  spin = particle.getSpinLikeQuantumNumber(QuantumNumberIDs::SPIN);
  spin_sum -= spin.J_numerator_ / spin.J_denominator_;

  return std::pow(-1.0, spin_sum) * parity_product;
}

void TopologyAmplitudeFactory::appendParticleInfoToName(std::stringstream &name,
    const ParticleStateInfo &state_info) const {
  name << state_info.pid_information_.name_ << "_J";

  name
      << state_info.spin_information_.J_numerator_
          / state_info.spin_information_.J_denominator_ << "_lambda";
  name
      << 1.0 * state_info.spin_information_.J_z_numerator_
          / state_info.spin_information_.J_denominator_;
}

std::pair<std::shared_ptr<DoubleParameter>, bool> TopologyAmplitudeFactory::generateGlobalParameter(
    const boost::property_tree::ptree& pt, const std::string& name,
    const std::string& related_name) {
  auto temp_param = generateDoubleParameter(pt, name);

  // check list for this parameter name
  auto result = global_parameter_list_.find(name);
  bool used_related_parameter(false);
  bool create_new(false);

  if (result == global_parameter_list_.end()) {
    if (parity_conserved && related_name != "") {
      // check parameter map for parity related phase
      auto result = global_parameter_list_.find(related_name);
      if (result != global_parameter_list_.end())
        used_related_parameter = true;
      else
        create_new = true;
    }
    else
      create_new = true;
  }

  if (create_new) {
    result =
        global_parameter_list_.insert(std::make_pair(name, temp_param)).first;
  }
  else if (used_related_parameter) {
    // ok check if both parameters have same fixing value
    if (result->second->IsFixed() != temp_param->IsFixed()) {
      if (temp_param->IsFixed()) {
        result->second->FixParameter(true);
      }
      BOOST_LOG_TRIVIAL(warning)<<"TopologyAmplitudeFactory::generateGlobalParameter(): "
      <<"two identical two body decays do not carry the same strength and phase fix property";
    }
  }

  return std::make_pair(result->second, used_related_parameter);
}

void TopologyAmplitudeFactory::fixParticlePropertyTree(
    boost::property_tree::ptree& pt, const IDInfo& id_info) const {
  // first create name for this id_info
  if (!pt.get_child_optional("mass").is_initialized()) {
    auto particle_properties = ComPWA::PhysConst::Instance().findParticle(
        id_info.particle_id_);
    pt.put("mass.value", particle_properties.mass_);
    pt.put("mass.fix", 1);
    pt.put("mass.min", 0.5 * particle_properties.mass_);
    pt.put("mass.max", 1.5 * particle_properties.mass_);
    pt.put("width.value", particle_properties.width_);
    pt.put("width.fix", 1);
    pt.put("width.min", 0.5 * particle_properties.width_);
    pt.put("width.max", 10.0 * particle_properties.width_);
  }
}

std::shared_ptr<DoubleParameter> TopologyAmplitudeFactory::generateDoubleParameter(
    const boost::property_tree::ptree& pt, const std::string& name) const {
  std::shared_ptr<DoubleParameter> new_parameter(
      new DoubleParameter(name, pt.get<double>("value")));
  if (pt.get<double>("min") != 0.0 || pt.get<double>("max") != 0.0)
    new_parameter->SetMinMax(pt.get<double>("min"), pt.get<double>("max"));
  new_parameter->FixParameter(pt.get<bool>("fix"));

  return new_parameter;
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
} /* namespace Physics */
} /* namespace ComPWA */
