// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

// Define Boost test module
#define BOOST_TEST_MODULE HelicityFormalism

#include "Core/Logging.hpp"
#include "Core/FourMomentum.hpp"
#include "Core/Properties.hpp"
#include "Core/Utils.hpp"
#include "Data/DataSet.hpp"
#include "Data/Generate.hpp"
#include "Data/Root/RootGenerator.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <vector>

using ComPWA::Physics::SubSystem;

std::ostream &operator<<(std::ostream &os, const std::vector<double> &v) {
  for (auto x : v) {
    os << x << ", ";
  }
  return os;
}

/// Calculation of helicity angle.
/// See (Martin and Spearman, Elementary Particle Theory. 1970)
double helicityAngle(double M, double m, double m2, double mSpec,
                     double invMassSqA, double invMassSqB) {
  // Calculate energy and momentum of m1/m2 in the invMassSqA rest frame
  double CmdEnergy =
      (invMassSqA + m * m - m2 * m2) / (2.0 * std::sqrt(invMassSqA));
  double CmsMomentum = CmdEnergy * CmdEnergy - m * m;
  // Calculate energy and momentum of mSpec in the invMassSqA rest frame
  double eSpecCms =
      (M * M - mSpec * mSpec - invMassSqA) / (2.0 * std::sqrt(invMassSqA));
  double pSpecCms = eSpecCms * eSpecCms - mSpec * mSpec;
  double CosAngle =
      -(invMassSqB - m * m - mSpec * mSpec - 2.0 * CmdEnergy * eSpecCms) /
      (2.0 * std::sqrt(CmsMomentum * pSpecCms));

  return CosAngle;
}

// Define Boost test suite (no idea what's the difference to TEST_MODULE)
BOOST_AUTO_TEST_SUITE(HelicityFormalism)

const std::string HelicityTestParticles = R"####(
<ParticleList>
  <Particle Name='K+'>
    <Pid>321</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_chargedKaon'>
      <Value>0.493677</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='1'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Spin' Type='IsoSpin' Value='0.5' Projection='0.5'/>
  </Particle>
  <Particle Name='K-'>
    <Pid>-321</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_chargedKaon'>
      <Value>0.493677</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='-1'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Spin' Type='IsoSpin' Value='0.5' Projection='-0.5'/>
  </Particle>
  <Particle Name='K_S0'>
    <Pid>310</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_neutralKaon'>
      <Value>0.497614</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Spin' Type='IsoSpin' Value='0.5' Projection='0.5'/>
  </Particle>
  <Particle Name='D0'>
    <Pid>421</Pid>
    <Parameter Type='Mass' Name='Mass_D0'>
      <Value>1.86484</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <DecayInfo Type='relativisticBreitWigner'>
      <FormFactor Type='0' />
      <Parameter Class='Double' Type='Mass' Name='Width_D0'>
        <Value>0.000623</Value>
        <Fix>true</Fix>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_D0'>
        <Value>2.5</Value>
        <Fix>true</Fix>
        <Min>2.0</Min>
        <Max>3.0</Max>
      </Parameter>
    </DecayInfo>
  </Particle>
</ParticleList>
)####";

/*! Test application for the calculation of the helicity angle.
 * As example the decay D0->KsK-K+ is used.
 * Some test may be specific to this decay but if all test are passed we can
 * be sure that helicity angles are calculated correctly.
 *
 * A note on numerical precision:
 * We compare all valued using float precision only since otherwise variables
 * calculated using different methods are not completely equal on double
 * precision. Every not an then even with float precision the test fails.
 */
BOOST_AUTO_TEST_CASE(HelicityAngleTest) {
  ComPWA::Logging Log("", "debug");

  std::stringstream ModelStream;
  // Construct particle list from XML tree
  ModelStream << HelicityTestParticles;
  auto partL = ComPWA::readParticles(ModelStream);

  // Construct HelicityKinematics by hand
  std::vector<ComPWA::pid> FinalState, InitialState;
  InitialState.push_back(421);
  FinalState.push_back(310);
  FinalState.push_back(-321);
  FinalState.push_back(321);
  auto kin =
      std::make_shared<ComPWA::Physics::HelicityFormalism::HelicityKinematics>(
          partL, InitialState, FinalState);

  // Generate phsp sample
  ComPWA::Data::Root::RootGenerator Generator(
      kin->getParticleStateTransitionKinematicsInfo());

  ComPWA::Data::Root::RootUniformRealGenerator RandomGenerator(123);

  auto Sample = ComPWA::Data::generatePhsp(20, Generator, RandomGenerator);

  bool useDerivedMassSq = false;

  ComPWA::Event ev;
  // We add an event of D0->KsK-K+ from data. Since the decay KS->pipi is not
  // constraint to the KS mass we use the relation:
  // (sqrtS * sqrtS + m1 * m1 + m2 * m2 + m3 * m3 - m23sq - m13sq)
  // to calculated the third invariant mass. The results for the sub systems
  // listed below are:
  // m23sq=1.4014 m13sq=1.52861 m12sq=1.28267
  // cosTheta12_23=0.151776 cosTheta12_CP=-0.178456
  // cosTheta13_12=0.178456 cosTheta13_CP=-0.151776
  // cosTheta23_12=0.318779 cosTheta23_CP=-0.318779
  // Angles are symmetric! Set @useDerivedMassSq = true.
  //
  // In case all invariant masses are calculated independently
  // (@useDerivedMassSq = false)the angles are not symmetric anymore.
  // In this case a couple of tests are supposed to fail.
  // m23sq=1.4014 m13sq=1.52861 m12sq=1.26798
  // cosTheta12_23=0.171554 cosTheta12_CP=-0.178456
  // cosTheta13_12=0.219724 cosTheta13_CP=-0.13337
  // cosTheta23_12=0.356843 cosTheta23_CP=-0.318779
  //    ev.addParticle(ComPWA::Particle(
  //        std::array<double, 4>{{-0.00827061, -0.242581, -0.335833,
  //        0.636104}},
  //        310));
  //    ev.addParticle(ComPWA::Particle(
  //        std::array<double, 4>{{-0.158637, -0.149132, 0.199913, 0.575405}},
  //        321));
  //    ev.addParticle(ComPWA::Particle(
  //        std::array<double, 4>{{-0.0236227, 0.453598, -0.0330521, 0.671656}},
  //        -321));
  //    sample->pushEvent(ev);

  double m1 = findParticle(partL, FinalState.at(0)).getMass().Value;
  double m2 = findParticle(partL, FinalState.at(1)).getMass().Value;
  double m3 = findParticle(partL, FinalState.at(2)).getMass().Value;
  double sqrtS = findParticle(partL, InitialState.at(0)).getMass().Value;

  // Define SubSystems that we want to test. The systems denoted my *_CP are
  // the CP conjugate systems the are used in the relation between D0 amplitude
  // A and D0bar amplitude Abar:
  // A(m_12^2,m_13^2) = Abar(m_13^2, m_12^2) -> A(sys12) = Abar(sys12_CP)
  // This is very specific to this decay.

  SubSystem sys12({{0}, {1}}, {2}, {});
  kin->registerSubSystem(sys12);
  SubSystem sys12_CP({{0}, {2}}, {1}, {});
  kin->registerSubSystem(sys12_CP);

  SubSystem sys13({{2}, {0}}, {1}, {});
  kin->registerSubSystem(sys13);
  SubSystem sys13_CP({{1}, {0}}, {2}, {});
  kin->registerSubSystem(sys13_CP);

  SubSystem sys23({{2}, {1}}, {0}, {});
  kin->registerSubSystem(sys23);
  SubSystem sys23_CP({{1}, {2}}, {0}, {});
  kin->registerSubSystem(sys23_CP);

  LOG(INFO) << "Loop over phsp events....";
  for (auto Event : Sample.Events) {
    double m23sq =
        (Event.FourMomenta[1] + Event.FourMomenta[2]).invariantMassSquared();
    double m13sq =
        (Event.FourMomenta[0] + Event.FourMomenta[2]).invariantMassSquared();
    double m12sq;
    if (useDerivedMassSq)
      m12sq = (sqrtS * sqrtS + m1 * m1 + m2 * m2 + m3 * m3 - m23sq - m13sq);
    else
      m12sq =
          (Event.FourMomenta[0] + Event.FourMomenta[1]).invariantMassSquared();

    double RelativeTolerance(1e-6);
    //------------ Restframe (12) -------------
    // Angle in the rest frame of (12) between (1) and (3)
    double cosTheta12_13 = helicityAngle(sqrtS, m1, m2, m3, m12sq, m13sq);
    double cosTheta12_13_2 = helicityAngle(sqrtS, m2, m1, m3, m12sq, m23sq);
    // Equal to the same angle calculated from m23sq
    BOOST_CHECK_CLOSE(cosTheta12_13, -1 * cosTheta12_13_2, RelativeTolerance);

    // Angle in the rest frame of (12) between (2) and (3)
    double cosTheta12_23 = helicityAngle(sqrtS, m2, m1, m3, m12sq, m23sq);
    double cosTheta12_23_2 = helicityAngle(sqrtS, m1, m2, m3, m12sq, m13sq);
    // Equal to the same angle calculated from m13sq
    BOOST_CHECK_CLOSE(cosTheta12_23, -1 * cosTheta12_23_2, RelativeTolerance);

    BOOST_CHECK_CLOSE(cosTheta12_13, -1 * cosTheta12_23, RelativeTolerance);
    BOOST_CHECK_CLOSE(cosTheta12_13, -1 * cosTheta12_23, RelativeTolerance);
    BOOST_CHECK_CLOSE(cosTheta12_13_2, -1 * cosTheta12_23_2, RelativeTolerance);

    //------------ Restframe (13) -------------
    // Angle in the rest frame of (13) between (1) and (2)
    double cosTheta13_12 = helicityAngle(sqrtS, m1, m3, m2, m13sq, m12sq);
    double cosTheta13_12_2 = helicityAngle(sqrtS, m3, m1, m2, m13sq, m23sq);
    BOOST_CHECK_CLOSE(cosTheta13_12, -1 * cosTheta13_12_2, RelativeTolerance);
    // Angle in the rest frame of (13) between (3) and (2)
    double cosTheta13_23 = helicityAngle(sqrtS, m3, m1, m2, m13sq, m23sq);
    double cosTheta13_23_2 = helicityAngle(sqrtS, m1, m3, m2, m13sq, m12sq);
    BOOST_CHECK_CLOSE(cosTheta13_23, -1 * cosTheta13_23_2, RelativeTolerance);
    BOOST_CHECK_CLOSE(cosTheta13_12, -1 * cosTheta13_23, RelativeTolerance);
    BOOST_CHECK_CLOSE(cosTheta13_12_2, -1 * cosTheta13_23_2, RelativeTolerance);

    //------------ Restframe (23) -------------
    // Angle in the rest frame of (23) between (2) and (1)
    double cosTheta23_12 = helicityAngle(sqrtS, m2, m3, m1, m23sq, m12sq);
    double cosTheta23_12_2 = helicityAngle(sqrtS, m3, m2, m1, m23sq, m13sq);
    BOOST_CHECK_CLOSE(cosTheta23_12, -1 * cosTheta23_12_2, RelativeTolerance);

    // Angle in the rest frame of (23) between (3) and (1)
    double cosTheta23_13 = helicityAngle(sqrtS, m3, m2, m1, m23sq, m13sq);
    double cosTheta23_13_2 = helicityAngle(sqrtS, m2, m3, m1, m23sq, m12sq);
    BOOST_CHECK_CLOSE(cosTheta23_13, -1 * cosTheta23_13_2, RelativeTolerance);

    BOOST_CHECK_CLOSE(cosTheta23_12, -1 * cosTheta23_13, RelativeTolerance);

    //------------- Test of HelicityKinematics -----------------
    // Check if calculation of helicity angles corresponds to the previously
    // calculated values

    LOG(DEBUG) << "-------- NEW EVENT ----------";
    auto p12 = kin->calculateHelicityAngles(Event, sys12);
    auto p12_CP = kin->calculateHelicityAngles(Event, sys12_CP);
    // BOOST_CHECK_EQUAL((float)p12.value(1), (float)cosTheta12_23);
    BOOST_CHECK(ComPWA::Utils::equal(std::cos(p12.first), cosTheta12_23, 1000));

    LOG(DEBUG) << "-------- (12) ----------";
    LOG(DEBUG) << sys12 << " : " << p12.first << ", " << p12.second;
    LOG(DEBUG) << sys12_CP << " : " << p12_CP.first << ", " << p12_CP.second;
    LOG(DEBUG) << "cosTheta12 "
               << helicityAngle(sqrtS, m2, m1, m3, m12sq, m23sq)
               << " CP: " << helicityAngle(sqrtS, m3, m1, m2, m13sq, m23sq);

    auto p13 = kin->calculateHelicityAngles(Event, sys13);
    auto p13_CP = kin->calculateHelicityAngles(Event, sys13_CP);
    BOOST_CHECK_CLOSE(std::cos(p13.first), cosTheta13_12, RelativeTolerance);

    BOOST_CHECK_CLOSE(std::cos(p13.first), -1 * std::cos(p12_CP.first),
                      RelativeTolerance);
    BOOST_CHECK_CLOSE(std::cos(p12.first), -1 * std::cos(p13_CP.first),
                      RelativeTolerance);

    LOG(DEBUG) << "-------- (13) ----------";
    LOG(DEBUG) << sys13 << " : " << p13.first << ", " << p13.second;
    LOG(DEBUG) << sys13_CP << " : " << p13_CP.first << ", " << p13_CP.second;
    LOG(DEBUG) << "cosTheta13 "
               << helicityAngle(sqrtS, m1, m3, m2, m13sq, m12sq)
               << " CP: " << helicityAngle(sqrtS, m1, m2, m3, m12sq, m13sq);

    auto p23 = kin->calculateHelicityAngles(Event, sys23);
    auto p23_CP = kin->calculateHelicityAngles(Event, sys23_CP);
    BOOST_CHECK_CLOSE(std::cos(p23.first), cosTheta23_12, RelativeTolerance);
    BOOST_CHECK_CLOSE(std::cos(p23.first), -1 * std::cos(p23_CP.first),
                      RelativeTolerance);

    LOG(DEBUG) << "-------- (23) ----------";
    LOG(DEBUG) << sys23 << " : " << p23.first << ", " << p23.second;
    LOG(DEBUG) << sys23_CP << " : " << p23_CP.first << ", " << p23_CP.second;
    LOG(DEBUG) << "cosTheta23 "
               << helicityAngle(sqrtS, m2, m3, m1, m23sq, m12sq)
               << " CP: " << helicityAngle(sqrtS, m3, m2, m1, m23sq, m13sq);
  }
}

BOOST_AUTO_TEST_SUITE_END()
