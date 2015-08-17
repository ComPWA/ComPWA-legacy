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

#include "HelicityKinematics.hpp"
#include "FinalStateParticleCombinatorics.hpp"
#include "Core/Event.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Particle.hpp"

namespace HelicityFormalism {

HelicityKinematics::HelicityKinematics() {
}

HelicityKinematics::~HelicityKinematics() {
}

bool HelicityKinematics::isWithinPhsp(const dataPoint& point) {
}

double HelicityKinematics::getMotherMass() const {
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

   gsl_monte_function F = { &phspFunc, dim, const_cast<DalitzKinematics*>(this) };

   gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(dim);
   gsl_monte_vegas_integrate(&F, xLimit_low, xLimit_high, 2, calls, r, s, &res,
   &err);
   gsl_monte_vegas_free(s);
   BOOST_LOG_TRIVIAL(debug)
   << "HelicityKinematics::calculatePSArea() phase space area (MC integration): " << "("
   << res << "+-" << err << ") relAcc [%]: " << 100 * err / res;

   PS_area_ = res;
   is_PS_area_calculated_ = 1;
   return;*/
}

void HelicityKinematics::setDecayTopologies(
    const std::vector<DecayTopology>& decay_topologies) {
  decay_topologies_ = decay_topologies;
}

void HelicityKinematics::init(const Event& event) {
  FinalStateParticleCombinatorics fsp_combinatorics;
  fsp_combinatorics.init(createFSParticleList(), event);

  // now for each decay topology  create a index list
  std::vector<DecayTopology>::const_iterator decay_topology_iter;
  for (decay_topology_iter = decay_topologies_.begin();
      decay_topology_iter != decay_topologies_.end(); ++decay_topology_iter) {

    std::vector<IndexMapping> mappings =
        fsp_combinatorics.getUniqueParticleMappingsSubsetForTopology(
            *decay_topology_iter);

    std::vector<IndexList> data_point_index_list;
    for (unsigned int mapping_index = 0; mapping_index < mappings.size();
        ++mapping_index) {
      data_point_index_list.push_back(
          createDataPointIndexListForTopology(*decay_topology_iter,
              mappings[mapping_index]));
    }

    topology_amplitude_data_point_index_lists_.push_back(data_point_index_list);
  }
}

std::vector<ParticleStateInfo> HelicityKinematics::createFSParticleList() const {
  std::vector<ParticleStateInfo> final_state_particle_pool;
  // just take the first topology (all should have the same fs particle content)
  if (decay_topologies_.size() > 0) {
    const DecayTopology& decay_topology = decay_topologies_.front();

    std::set<ParticleStateInfo> temp_particle_set;

    for (auto final_state_particle_lists_iter =
        decay_topology.final_state_content_lists_.begin();
        final_state_particle_lists_iter
            != decay_topology.final_state_content_lists_.end();
        ++final_state_particle_lists_iter) {
      temp_particle_set.insert(final_state_particle_lists_iter->begin(),
          final_state_particle_lists_iter->end());
    }
  }
  return final_state_particle_pool;
}

IndexList HelicityKinematics::createDataPointIndexListForTopology(
    const DecayTopology& topology,
    const IndexMapping& fs_particle_mapping) {

  IndexList topology_amplitude_data_point_index_list;
  // loop through the decay topology
  for (auto decay_node = topology.decay_node_infos_.begin();
      decay_node != topology.decay_node_infos_.end(); ++decay_node) {

    TwoBodyDecayIndices decay_indices;
    // for each decay node: create and index list of final state particles
    // for all occuring cms frames (maybe they already exist... so check)
    if (0 < decay_node->mother_index_) {
      decay_indices.mother_index_ = convertAndStoreParticleList(
          topology.final_state_content_lists_[decay_node->mother_index_],
          fs_particle_mapping);
    }
    else {
      decay_indices.mother_index_ = decay_node->mother_index_;
    }
    decay_indices.decay_products_.first = convertAndStoreParticleList(
        topology.final_state_content_lists_[decay_node->decay_products_.first],
        fs_particle_mapping);
    decay_indices.decay_products_.second = convertAndStoreParticleList(
        topology.final_state_content_lists_[decay_node->decay_products_.second],
        fs_particle_mapping);

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
    topology_amplitude_data_point_index_list.push_back(index_pair_index);
  }
  return topology_amplitude_data_point_index_list;
}

unsigned int HelicityKinematics::convertAndStoreParticleList(
    const std::vector<ParticleStateInfo>& particle_list,
    const IndexMapping& fs_particle_mapping) {

  IndexList event_particle_index_list = convertParticleListToEventIndexList(
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

IndexList HelicityKinematics::convertParticleListToEventIndexList(
    const std::vector<ParticleStateInfo>& particle_list,
    const IndexMapping& fs_particle_mapping) const {
  IndexList event_particle_index_list;
  event_particle_index_list.reserve(particle_list.size());
  for (auto particle_iter = particle_list.begin();
      particle_iter != particle_list.end(); ++particle_iter) {
    event_particle_index_list.push_back(
        fs_particle_mapping.find(particle_iter->id_)->second);
  }
  return event_particle_index_list;
}

void HelicityKinematics::translateEventToDataPoint(const Event& event,
    dataPoint& point) const {
  point.reset(3 * unique_occurring_decay_product_index_list_links_.size());
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
}

// IMPORTANT: the ordering of the 4 vectors is exactly the same as in the
// unique_occurring_decay_event_fs_index_lists_
std::vector<Vector4<double> > HelicityKinematics::createRequired4Vectors(
    const Event& event) const {
  std::vector<Vector4<double> > unique_occurring_cms_4vectors;
  unique_occurring_cms_4vectors.reserve(
      unique_occurring_decay_event_fs_index_lists_.size());

  for (unsigned int fs_particle_lists_index = 0;
      fs_particle_lists_index
          < unique_occurring_decay_event_fs_index_lists_.size();
      ++fs_particle_lists_index) {
    Vector4<double> temp_vector;
    for (unsigned int event_particle_list_index = 0;
        event_particle_list_index
            < unique_occurring_decay_event_fs_index_lists_[fs_particle_lists_index].size();
        ++event_particle_list_index) {
      addParticleToCMS4Vector(
          event.getParticle(
              unique_occurring_decay_event_fs_index_lists_[fs_particle_lists_index][event_particle_list_index]),
          temp_vector);
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
  Vector4<double> particle1_4vector(
      unique_occurring_cms_4vectors[two_body_state_indices.decay_products_.first]);
  Vector4<double> particle2_4vector(
      unique_occurring_cms_4vectors[two_body_state_indices.decay_products_.second]);
  // define the two body state
  Vector4<double> decaying_state(
      unique_occurring_cms_4vectors[two_body_state_indices.decay_products_.first]
          + unique_occurring_cms_4vectors[two_body_state_indices.decay_products_.second]);

  // at first add the invariant mass squared to the data point
  point.setVal(data_point_fill_position++, decaying_state.Mass2());
  // then add the invariant masses of the daughters
  point.setVal(data_point_fill_position++, particle1_4vector.Mass());
  point.setVal(data_point_fill_position++, particle2_4vector.Mass());

  // boost particle1 into the rest frame of the two body state
  particle1_4vector.Boost(decaying_state);

  // if this decay node is the not the top level decay we have to do some
  // boosting
  if (two_body_state_indices.mother_index_ > 0) {
    // define mother state
    Vector4<double> mother(
        unique_occurring_cms_4vectors[two_body_state_indices.mother_index_]);
    // then boost the two body state into the rest frame of its mother
    decaying_state.Boost(mother);
    // now determine the theta and phi values of the boosted particle1 vector
    // with respect to the boosted two body state
    double rotation_phi = decaying_state.Phi();
    particle1_4vector.Rotate(rotation_phi, decaying_state.Theta(),
        -rotation_phi);
  }
  // now just get the theta and phi angles of the boosted particle 1
  point.setVal(data_point_fill_position++, particle1_4vector.Theta());
  point.setVal(data_point_fill_position++, particle1_4vector.Phi());
}

std::vector<std::vector<IndexList> > HelicityKinematics::getTopologyAmplitudeDataPointIndexLists() const {
  return topology_amplitude_data_point_index_lists_;
}

double HelicityKinematics::getMass(unsigned int num) const {
}

double HelicityKinematics::getMass(std::string name) const {
}

} /* namespace HelicityFormalism */
