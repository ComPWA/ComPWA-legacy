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

// this include is required to define the size_t type
// needed by the gsl monte carlo stuff
#include <cstddef>

#include "gsl/gsl_monte.h"
#include "gsl/gsl_monte_vegas.h"

#include "Core/Event.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Particle.hpp"
#include "Core/PhysConst.hpp"

#include "Physics/HelicityAmplitude/HelicityKinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

HelicityKinematics::HelicityKinematics() {
}

HelicityKinematics::~HelicityKinematics() {
}

bool HelicityKinematics::isWithinPhsp(const dataPoint& point) {
  for (unsigned int i = 0;
      i < unique_occurring_decay_product_index_list_links_.size(); ++i) {
    double mother_mass(sqrt(point.getVal((i * getNumberOfVariables()))));
    double particle1_mass(point.getVal((i * getNumberOfVariables()) + 1));
    double particle2_mass(point.getVal((i * getNumberOfVariables()) + 2));

    if (mother_mass < particle1_mass + particle2_mass)
      return false;
  }

  return true;
}

double HelicityKinematics::getMotherMass() const {
  return mother_mass_;
}

double HelicityKinematics::calculatePSArea() {
  /*size_t dim = 2;
   double res = 0.0, err = 0.0;

   //set limits: we assume that x[0]=m13sq and x[1]=m23sq
   double xLimit_low[2] = { m13_sq_min, m23_sq_min };
   double xLimit_high[2] = { m13_sq_max, m23_sq_max };

   size_t calls = 2000000;
   gsl_rng_env_setup();
   const gsl_rng_type *T = gsl_rng_default;    //type of random generator
   gsl_rng *r = gsl_rng_alloc(T);    //random generator

   gsl_monte_function F = { &heliphspFunc, dim, const_cast<HelicityKinematics*>(this) };

   gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(dim);
   gsl_monte_vegas_integrate(&F, xLimit_low, xLimit_high, 2, calls, r, s, &res,
   &err);
   gsl_monte_vegas_free(s);

   BOOST_LOG_TRIVIAL(debug)
   << "HelicityKinematics::calculatePSArea() phase space area (MC integration): " << "("
   << res << "+-" << err << ") relAcc [%]: " << 100 * err / res;*/

  /*TGenPhaseSpace gen;
   TLorentzVector mother;
   double *masses = {};
   gen.SetDecay(mother, num_particles, masses);

   double precision(1e-5);
   unsigned int max_calls(10000000);

   double added_weights(0.0);
   unsigned int counter(0);
   double integral(1.0);
   double previous_value(0.0);

   while(std::fabs((integral-previous_value)/integral) > precision && counter < max_calls) {
   previous_value = integral;
   added_weights += gen.Generate();
   ++counter;
   integral = added_weights/counter;
   }

   return integral;*/
  return 1.0;
}

void HelicityKinematics::setDecayTopologies(
    const std::vector<TwoBodyDecayTopology>& decay_topologies) {
  decay_topologies_ = decay_topologies;

  if (decay_topologies_.size() > 0) {
    const TwoBodyDecayTopology& dt = decay_topologies_[0];

    const std::vector<IDInfo> &fs_particle_list = dt.final_state_id_list_;

    ComPWA::PhysConst& physics_constants = ComPWA::PhysConst::Instance();

    mother_mass_ =
        physics_constants.findParticle(dt.top_node_id_info_.name_).mass_;

    number_of_particles_ = fs_particle_list.size();
    // now set the particle mass and name vectors
    for (auto fs_particle = fs_particle_list.begin();
        fs_particle != fs_particle_list.end(); ++fs_particle) {
      particle_names_.push_back(fs_particle->name_);
      masses_.push_back(
          physics_constants.findParticle(fs_particle->name_).mass_);
    }
  }
}

void HelicityKinematics::init(
    const FinalStateParticleCombinatorics& fsp_combinatorics) {

  // now add the variable names needed for one decay vertex
  variable_names_.push_back("cms_mass_squared");
  variable_names_.push_back("daughter1_mass");
  variable_names_.push_back("daughter2_mass");
  variable_names_.push_back("helicity_angle_theta");
  variable_names_.push_back("helicity_angle_phi");

  fsp_combinatorics_ = fsp_combinatorics;
  // now for each decay topology create a index list
  for (auto const& decay_topology : decay_topologies_) {

    std::vector<IndexMapping> mappings =
        fsp_combinatorics_.getUniqueParticleMappingsSubsetForTopology(
            decay_topology);

    std::vector<IndexList> data_point_index_list;
    for (auto const& mapping : mappings) {

      IndexList topology_amplitude_data_point_index_list;

      for (auto evalution_order_index : decay_topology.unique_id_decay_node_order_) {
        buildDataPointIndexListForTopology(evalution_order_index,
            decay_topology, mapping, topology_amplitude_data_point_index_list);
      }

      data_point_index_list.push_back(topology_amplitude_data_point_index_list);
    }

    topology_amplitude_data_point_index_lists_.push_back(data_point_index_list);
  }
}

void HelicityKinematics::buildDataPointIndexListForTopology(unsigned int index,
    const TwoBodyDecayTopology& topology,
    const IndexMapping& fs_particle_mapping, IndexList& data_point_index_list) {

  // if its a leaf then just return
  if (topology.particle_unique_id_decay_tree_.find(index)
      == topology.particle_unique_id_decay_tree_.end())
    return;

  const IndexList& daughter_unique_id_list =
      topology.particle_unique_id_decay_tree_.at(index);

  TwoBodyDecayIndices decay_indices;

  decay_indices.decay_state_index_ = convertAndStoreParticleIndexList(
      topology.final_state_content_unique_id_mapping_.at(index),
      fs_particle_mapping);

  unsigned int mother_index(topology.findMotherIndex(index));
  if (mother_index != index) {
    decay_indices.mother_index_ = convertAndStoreParticleIndexList(
        topology.final_state_content_unique_id_mapping_.at(mother_index),
        fs_particle_mapping);
  }
  else {
    decay_indices.mother_index_ = decay_indices.decay_state_index_;
  }

  unsigned int final_state_counter(0);
  for (auto daughter_unique_id = daughter_unique_id_list.begin();
      daughter_unique_id != daughter_unique_id_list.end();
      ++daughter_unique_id) {

    if (final_state_counter == 0)
      decay_indices.decay_products_.first = convertAndStoreParticleIndexList(
          topology.final_state_content_unique_id_mapping_.at(
              *daughter_unique_id), fs_particle_mapping);
    else
      decay_indices.decay_products_.second = convertAndStoreParticleIndexList(
          topology.final_state_content_unique_id_mapping_.at(
              *daughter_unique_id), fs_particle_mapping);

    ++final_state_counter;
  }

  // then also link these cms frames to the decay product index links
  // (again the may already exist)
  auto found_iter = std::find(
      unique_occurring_decay_product_index_list_links_.begin(),
      unique_occurring_decay_product_index_list_links_.end(), decay_indices);

  unsigned int index_pair_index;
  if (found_iter != unique_occurring_decay_product_index_list_links_.end()) {
    index_pair_index = found_iter
        - unique_occurring_decay_product_index_list_links_.begin();
  }
  else {
    unique_occurring_decay_product_index_list_links_.push_back(decay_indices);
    index_pair_index = unique_occurring_decay_product_index_list_links_.size()
        - 1;
  }

  // finally add the index, which points to this
  // decay fs event index list pair, to the index list
  data_point_index_list.push_back(getNumberOfVariables() * index_pair_index);

  // and recurse to daughters
  /*for (auto daughter_unique_id = daughter_unique_id_list.begin();
   daughter_unique_id != daughter_unique_id_list.end();
   ++daughter_unique_id) {
   recursivelyBuildDataPointIndexListForTopology(*daughter_unique_id, topology,
   fs_particle_mapping, data_point_index_list);
   }*/
}

unsigned int HelicityKinematics::convertAndStoreParticleIndexList(
    const IndexList& particle_list, const IndexMapping& fs_particle_mapping) {

  IndexList event_particle_index_list = convertParticleIDToEventIndexList(
      particle_list, fs_particle_mapping);

  auto found_iter = std::find(
      unique_occurring_decay_event_fs_index_lists_.begin(),
      unique_occurring_decay_event_fs_index_lists_.end(),
      event_particle_index_list);

  if (found_iter != unique_occurring_decay_event_fs_index_lists_.end()) {
    return found_iter - unique_occurring_decay_event_fs_index_lists_.begin();
  }
  else {
    unique_occurring_decay_event_fs_index_lists_.push_back(
        event_particle_index_list);
    return unique_occurring_decay_event_fs_index_lists_.size() - 1;
  }
}

IndexList HelicityKinematics::convertParticleIDToEventIndexList(
    const IndexList& particle_id_list,
    const IndexMapping& fs_particle_mapping) const {
  IndexList event_particle_index_list;
  event_particle_index_list.reserve(particle_id_list.size());
  for (auto particle_iter = particle_id_list.begin();
      particle_iter != particle_id_list.end(); ++particle_iter) {
    event_particle_index_list.push_back(
        fs_particle_mapping.find(*particle_iter)->second);
  }

  std::sort(event_particle_index_list.begin(), event_particle_index_list.end());

  return event_particle_index_list;
}

void HelicityKinematics::translateEventToDataPoint(const Event& event,
    dataPoint& point) const {
  point.reset(
      getNumberOfVariables()
          * unique_occurring_decay_product_index_list_links_.size());
  // create all needed cms frame 4 vectors
  std::vector<Vector4<double> > unique_occurring_cms_4vectors =
      createRequired4Vectors(event);
  // then just loop over all occurring decays and determine the kinematic
  // variables and append them to the dataPoint in the same order
  unsigned int fill_position = 0;
  for (unsigned int decay_index = 0;
      decay_index < unique_occurring_decay_product_index_list_links_.size();
      ++decay_index) {
    fillPointWithBoostedKinematicVariables(point, unique_occurring_cms_4vectors,
        unique_occurring_decay_product_index_list_links_[decay_index],
        fill_position);
  }

  // check event to data point translation
  /*std::cout << "event:\n";
   for (unsigned int i = 0; i < event.getNParticles(); ++i) {
   std::cout << event.getParticle(i).E << ", " << event.getParticle(i).px
   << ", " << event.getParticle(i).py << ", " << event.getParticle(i).pz
   << std::endl;
   }
   std::cout << "4vectors:\n";
   for (unsigned int i = 0; i < unique_occurring_cms_4vectors.size(); ++i) {
   unique_occurring_cms_4vectors[i].Print(std::cout);
   }
   std::cout << "\ndata point:\n[";
   for (unsigned int i = 0; i < point.size(); ++i) {
   std::cout << point.getVal(i) << ", ";
   }
   std::cout << "]\n" << std::endl;*/
}

// IMPORTANT: the ordering of the 4 vectors is exactly the same as in the
// unique_occurring_decay_event_fs_index_lists_
std::vector<Vector4<double> > HelicityKinematics::createRequired4Vectors(
    const Event& event) const {
  std::vector<Vector4<double> > unique_occurring_cms_4vectors;
  // we create 3 vectors for each unique decay. Mother state + the two daughters = 3
  unique_occurring_cms_4vectors.reserve(
      3 * unique_occurring_decay_event_fs_index_lists_.size());

  for (auto const& fs_particle_list : unique_occurring_decay_event_fs_index_lists_) {
    Vector4<double> temp_vector(0.0, 0.0, 0.0, 0.0);

    if (fs_particle_list.size() == 1) {
      const Particle& p = event.getParticle(fs_particle_list[0]);
      temp_vector.SetP4(p.E, p.px, p.py, p.pz);
    }
    else {
      for (unsigned int event_particle_list_index = 0;
          event_particle_list_index < fs_particle_list.size();
          ++event_particle_list_index) {

        addParticleToCMS4Vector(
            event.getParticle(fs_particle_list[event_particle_list_index]),
            temp_vector);
      }
    }
    unique_occurring_cms_4vectors.push_back(temp_vector);
  }

  return unique_occurring_cms_4vectors;
}

void HelicityKinematics::addParticleToCMS4Vector(const Particle& event_particle,
    Vector4<double>& cms_4vector) const {
  double px(cms_4vector.Px());
  double py(cms_4vector.Py());
  double pz(cms_4vector.Pz());
  double E(cms_4vector.E());

  px += event_particle.px;
  py += event_particle.py;
  pz += event_particle.pz;
  E += event_particle.E;

  cms_4vector.SetP4(E, px, py, pz);
}

/**
 * Calculates theta, phi and the invariant mass square which are required for
 * the evaluation of the amplitudes in the helicity formalism.
 */
void HelicityKinematics::fillPointWithBoostedKinematicVariables(
    dataPoint& point,
    const std::vector<Vector4<double> >& unique_occurring_cms_4vectors,
    const TwoBodyDecayIndices& two_body_state_indices,
    unsigned int& data_point_fill_position) const {
  // define particle products of the two body decay
  Vector4<double> daughter1_4vector(
      unique_occurring_cms_4vectors[two_body_state_indices.decay_products_.first]);
  Vector4<double> daughter2_4vector(
      unique_occurring_cms_4vectors[two_body_state_indices.decay_products_.second]);
  // define the two body state
  Vector4<double> decaying_state(
      unique_occurring_cms_4vectors[two_body_state_indices.decay_state_index_]);

  // at first add the invariant mass squared to the data point
  point.setVal(data_point_fill_position, decaying_state.Mass2());
  // then add the invariant masses of the daughters
  point.setVal(++data_point_fill_position, daughter1_4vector.Mass());
  point.setVal(++data_point_fill_position, daughter2_4vector.Mass());

  // boost particle1 into the rest frame of the two body state
  daughter1_4vector.Boost(decaying_state);

  // if this decay node is the not the top level decay we have to do some
  // boosting
  if (two_body_state_indices.decay_state_index_
      != two_body_state_indices.mother_index_) {
    // define mother state
    Vector4<double> mother(
        unique_occurring_cms_4vectors[two_body_state_indices.mother_index_]);

    // then boost the two body state into the rest frame of its mother
    decaying_state.Boost(mother);
    // now determine the theta and phi values of the boosted particle1 vector
    // with respect to direction of flight of the boosted two body state
    // so we need to rotate
    //double rotation_theta = daughter1_4vector.Theta();
    //double rotation_phi = daughter1_4vector.Phi();
    //daughter1_4vector.RotateZ(-rotation_phi);
    //daughter1_4vector.RotateY(-rotation_theta);
    daughter1_4vector.Rotate(decaying_state.Phi(), decaying_state.Theta(), -decaying_state.Phi());
    //daughter1_4vector.RotateY(rotation_theta);
    //daughter1_4vector.RotateZ(rotation_phi);
  }

   /*Vector4<double> mother(std::sqrt(decaying_state.Mass2() + 1.0), 0., 0., 1.0);

   if (two_body_state_indices.decay_state_index_
   != two_body_state_indices.mother_index_) {
   mother =
   unique_occurring_cms_4vectors[two_body_state_indices.mother_index_];
   }

   decaying_state.Boost(mother);

   Vector4<double> decaying_state_rot(decaying_state);
   decaying_state_rot.RotateZ(-decaying_state.Phi());
   decaying_state_rot.RotateY(-decaying_state.Theta());

   particle1_4vector.Boost(mother);
   particle1_4vector.RotateZ(-decaying_state.Phi());
   particle1_4vector.RotateY(-decaying_state.Theta());

   particle1_4vector.Boost(decaying_state_rot);*/

  // now just get the theta and phi angles of the boosted particle 1
  point.setVal(++data_point_fill_position, daughter1_4vector.Theta());

  point.setVal(++data_point_fill_position, daughter1_4vector.Phi());
  ++data_point_fill_position;
}

std::vector<std::vector<IndexList> > HelicityKinematics::getTopologyAmplitudeDataPointIndexLists() const {
  return topology_amplitude_data_point_index_lists_;
}

double HelicityKinematics::getMass(unsigned int num) const {
  return masses_.at(num - 1);    // they start counting from 1, but i store from 0 on...
}

double HelicityKinematics::getMass(std::string name) const {
  unsigned int position = std::find(particle_names_.begin(),
      particle_names_.end(), name) - particle_names_.begin();
  return masses_.at(position);
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
