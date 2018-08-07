// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

// Define Boost test module
#define BOOST_TEST_MODULE HelicityFormalism

#include <vector>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>

#include "Core/Logging.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Particle.hpp"
#include "Core/Properties.hpp"
#include "DataReader/Data.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IncoherentIntensity.hpp"

#include "Tools/Generate.hpp"
#include "Tools/RootGenerator.hpp"

#include "qft++/Vector4.h"

using namespace ComPWA;
using namespace ComPWA::Physics::HelicityFormalism;

using ComPWA::QFT::Vector4;
using ComPWA::QFT::Tensor;

BOOST_AUTO_TEST_SUITE(HelicityFormalism)

std::pair<double, double>
pawianHelicityAngles(const Vector4<double> &motherRef,
                     const Vector4<double> &ref, const Vector4<double> &mother,
                     const Vector4<double> &daughter) {
  Vector4<double> result = daughter;
  Vector4<double> refTrafo = ref;
  Vector4<double> refRecoilTrafo = motherRef - ref;

  // boost all vectors into the mother rest frame
  result.Boost(mother);
  refTrafo.Boost(mother);
  refRecoilTrafo.Boost(mother);

  // rotate vectors so that refRecoilTrafo moves in the negative direction of
  // the z-axis
  result.RotateZ(-refTrafo.Phi());
  result.RotateY(TMath::Pi() - refTrafo.Theta());

  refRecoilTrafo.RotateZ(-refTrafo.Phi());
  refRecoilTrafo.RotateY(TMath::Pi() - refTrafo.Theta());

  // rotate around the z-axis so that refRecoil lies in the x-z plain
  result.RotateZ(TMath::Pi() - refRecoilTrafo.Phi());
  return std::make_pair(result.CosTheta(), result.Phi());
}

std::vector<std::pair<double, double>> allPawianHelicityAngles(
    Vector4<double> cms_state, Vector4<double> intermediate_state1,
    Vector4<double> intermediate_state2, Vector4<double> final) {
  std::vector<std::pair<double, double>> resultvec;

  Vector4<double> topref(0., 0., 0., -2.0);
  Vector4<double> top(std::sqrt(cms_state.Mass2() + 1.0), 0., 0., -1.0);

  resultvec.push_back(
      pawianHelicityAngles(topref, top, cms_state, intermediate_state1));
  resultvec.push_back(pawianHelicityAngles(top, cms_state, intermediate_state1,
                                           intermediate_state2));
  resultvec.push_back(pawianHelicityAngles(cms_state, intermediate_state1,
                                           intermediate_state2, final));

  return resultvec;
}

// ----- These functions have been copied from EvtGen and ported to qft++ -----
double evtGenDecayAngle(const Vector4<double> &p, const Vector4<double> &q,
                        const Vector4<double> &d) {

  double pd = p * d;
  double pq = p * q;
  double qd = q * d;
  double mp2 = p.Mass2();
  double mq2 = q.Mass2();
  double md2 = d.Mass2();

  double cost = (pd * mq2 - pq * qd) /
                sqrt((pq * pq - mq2 * mp2) * (qd * qd - mq2 * md2));

  return cost;
}

double mag2r3(const Vector4<double> &p0, const Vector4<double> &p1) {
  return std::pow(p0 * p1, 2) / p0.Mass2() - p1.Mass2();
}

// Calculate the 3-d dot product of 4-vectors p1 and p2 in the rest frame of
// 4-vector p0
double dotr3(const Vector4<double> &p0, const Vector4<double> &p1,
             const Vector4<double> &p2) {
  return 1 / p0.Mass2() * (p0 * p1) * (p0 * p2) - p1 * p2;
}

Tensor<double> dual(const Tensor<double> &t2) {
  Tensor<double> temp(2);
  temp(0, 0) = 0.0;
  temp(1, 1) = 0.0;
  temp(2, 2) = 0.0;
  temp(3, 3) = 0.0;

  temp(0, 1) = t2(3, 2) - t2(2, 3);
  temp(0, 2) = -t2(3, 1) + t2(1, 3);
  temp(0, 3) = t2(2, 1) - t2(1, 2);

  temp(1, 2) = -t2(3, 0) + t2(0, 3);
  temp(1, 3) = t2(2, 0) - t2(0, 2);

  temp(2, 3) = -t2(1, 0) + t2(0, 1);

  temp(1, 0) = -temp(0, 1);
  temp(2, 0) = -temp(0, 2);
  temp(3, 0) = -temp(0, 3);

  temp(2, 1) = -temp(1, 2);
  temp(3, 1) = -temp(1, 3);

  temp(3, 2) = -temp(2, 3);
  return temp;
}

// Calculate ( \vec{p1} cross \vec{p2} ) \cdot \vec{p3} in rest frame of object
double scalartripler3(const Vector4<double> &p0, const Vector4<double> &p1,
                      const Vector4<double> &p2, const Vector4<double> &p3) {
  Tensor<double> l = dual(p0 % p1) * p2;
  return -1.0 / p0.Mass() * (l * p3);
}

// Calculate phi using the given 4 vectors (all in the same frame)
double evtGenDecayAnglePhi(const Vector4<double> &z, const Vector4<double> &p,
                           const Vector4<double> &q, const Vector4<double> &d) {
  double eq = (p * q) / p.Mass();
  double ed = (p * d) / p.Mass();
  double mq = q.Mass();
  double q2 = mag2r3(p, q);
  double qd = dotr3(p, q, d);
  double zq = dotr3(p, z, q);
  double zd = dotr3(p, z, d);
  double alpha = (eq - mq) / (q2 * mq) * qd - ed / mq;

  double y = scalartripler3(p, z, q, d) + alpha * scalartripler3(p, z, q, q);
  double x = (zq * (qd + alpha * q2) - q2 * (zd + alpha * zq)) / std::sqrt(q2);

  double phi = std::atan2(y, x);

  return phi;
}

// ----- End EvtGen helper functions -----

std::pair<double, double> calculateEvtGenAngles(
    const Vector4<double> &grandparent, const Vector4<double> &parent,
    const Vector4<double> &decaying_state, const Vector4<double> &daughter) {
  return std::make_pair(
      evtGenDecayAngle(parent, decaying_state, daughter),
      evtGenDecayAnglePhi(grandparent, parent, decaying_state, daughter));
}

double calculatePhiDiff(double phi1, double phi2) {
  double phidiff = phi1 - phi2;
  if (phidiff > TMath::Pi())
    phidiff -= 2.0 * TMath::Pi();
  else if (phidiff < -TMath::Pi())
    phidiff += 2.0 * TMath::Pi();
  return phidiff;
}

/*!
 * This test case verifies the correctness of the theta and phi angles in a
 * complicated decay tree, by comparing the results with evt gen and the
 * pawian
 */
BOOST_AUTO_TEST_CASE(HelicityAnglesCorrectnessTest) {
  ComPWA::Logging log("", "debug");

  // Construct HelicityKinematics from XML tree
  boost::property_tree::ptree tr;
  boost::property_tree::xml_parser::read_xml(
      "HelicityAnglesCorrectnessTest-input.xml", tr);
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, tr);

  auto kin = std::make_shared<HelicityKinematics>(
      partL, tr.get_child("HelicityKinematics"));

  auto intensity = std::make_shared<ComPWA::Physics::IncoherentIntensity>(
      partL, kin, tr.get_child("Intensity"));

  // Generate phsp sample
  std::shared_ptr<ComPWA::Generator> gen(new ComPWA::Tools::RootGenerator(
      partL, kin->getKinematicsProperties().InitialState,
      kin->getKinematicsProperties().FinalState, 123));
  std::shared_ptr<ComPWA::DataReader::Data> sample(
      new ComPWA::DataReader::Data());

  ComPWA::Tools::generatePhsp(50, gen, sample);

  Vector4<double> top_vec4(0, 0, 0, 1);

  LOG(INFO) << "Loop over phsp events and comparison of angles....";
  for (auto ev : sample->events()) {
    DataPoint compwa_point;
    kin->convert(ev, compwa_point);

    // convert evt to evtgen 4 vectors
    std::vector<Vector4<double>> temp;
    for (auto const &part : ev.particles()) {
      temp.push_back(Vector4<double>(part.fourMomentum()));
    }

    Vector4<double> level0 = temp[0] + temp[1] + temp[2] + temp[3];
    Vector4<double> level1 = temp[1] + temp[2] + temp[3];
    Vector4<double> level2 = temp[2] + temp[3];
    Vector4<double> level3 = temp[2];

    // Pawian angles
    auto pawian_angles =
        allPawianHelicityAngles(level0, level1, level2, level3);

    // EvtGen angles
    std::vector<std::pair<double, double>> evtgen_angles;
    evtgen_angles.push_back(std::make_pair(level1.CosTheta(), level1.Phi()));
    evtgen_angles.push_back(
        calculateEvtGenAngles(top_vec4, level0, level1, level2));
    evtgen_angles.push_back(
        calculateEvtGenAngles(level0, level1, level2, level3));
    level1.Boost(level0);
    evtgen_angles[0] = std::make_pair(level1.CosTheta(), level1.Phi());

    // ComPWA angles
    std::vector<std::pair<double, double>> compwa_angles;
    compwa_angles.push_back(
        std::make_pair(compwa_point.values()[1], compwa_point.values()[2]));
    compwa_angles.push_back(
        std::make_pair(compwa_point.values()[4], compwa_point.values()[5]));
    compwa_angles.push_back(
        std::make_pair(compwa_point.values()[7], compwa_point.values()[8]));

    for (unsigned int i = 0; i < compwa_angles.size(); ++i) {
      BOOST_CHECK_EQUAL((float)compwa_angles[i].first,
                        (float)pawian_angles[i].first);
      BOOST_CHECK_EQUAL((float)compwa_angles[i].first,
                        (float)evtgen_angles[i].first);
      BOOST_CHECK_EQUAL((float)compwa_angles[i].second,
                        (float)evtgen_angles[i].second);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
