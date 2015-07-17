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

#include "gsl/gsl_monte.h"
#include "gsl/gsl_monte_vegas.h"

#include "HelicityKinematics.hpp"

namespace HelicityFormalism {

HelicityKinematics::HelicityKinematics(
    const std::vector<HelicityDecayTree>& decay_trees) :
    decay_trees_(decay_trees), is_PS_area_calculated_(false), PS_area_(0.0) {
}

HelicityKinematics::~HelicityKinematics() {
}

bool HelicityKinematics::isWithinPhsp(const dataPoint& point) const {
}

double HelicityKinematics::getMotherMass() const {
}

double HelicityKinematics::getPhspVolume() const {
  if (!is_PS_area_calculated_)
    calculatePSArea();
  return PS_area_;
}

void HelicityKinematics::calculatePSArea() {
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


void HelicityKinematics::generateFinalStateMapping(Event) {

}

void HelicityKinematics::eventToDataPoint(Event& ev, dataPoint& point) const {
  // we need some kind of event particle to helicity state particle matching
  // maybe a map of helicity decay particle ids to lists of final state particles
  // this part should be done only once in the beginning, since the event particle
  // ordering should be the same for each event.

  // loop over the different decay trees
  std::vector<HelicityDecayTree>::const_iterator decay_tree_iter;
  for (decay_tree_iter = decay_trees_.begin();
      decay_tree_iter != decay_trees_.end(); ++decay_tree_iter) {


    std::vector<ParticleStatePair> decay_products = decay_tree_iter->getNextDecayVertexConnectedFinalStateParticleList();

  }
}

/**
 * Calculates the boosted 4 vector of the first particle from the two body decay
 * in the helicity formalism.
 */
Vector4 HelicityKinematics::determineBoostedKinematicVariables(
    std::pair<Vector4, Vector4> two_body_state, Vector4 mother) {
  // define particle 1 of the two body decay
  PolVector particle1_4vector;
  particle1_4vector.SetP4(two_body_state.first);
  // define the two body state
  PolVector two_body_4vector;
  two_body_4vector.SetP4(two_body_state.first + two_body_state.second);
  // boost particle1 into the rest frame of the two body state
  particle1_4vector.Boost(two_body_rest_frame.GetP4());
  // then boost the two body state into the rest frame of its mother
  two_body_4vector.Boost(mother);
  // now determine the theta and phi values of the boosted particle1 vector
  // with respect to the boosted two body state

}

double HelicityKinematics::getMass(unsigned int num) const {
}

double HelicityKinematics::getMass(std::string name) const {
}

} /* namespace HelicityFormalism */
