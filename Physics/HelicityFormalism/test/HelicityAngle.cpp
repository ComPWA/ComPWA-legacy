// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

// Define Boost test module
#define BOOST_TEST_MODULE HelicityFormalism

#include <iostream>
#include <vector>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>

#include "Core/Logging.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Particle.hpp"
#include "Core/Properties.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

#include "Data/DataSet.hpp"
#include "Tools/Generate.hpp"
#include "Tools/RootGenerator.hpp"

using namespace ComPWA;
using namespace ComPWA::Physics::HelicityFormalism;
using ComPWA::Physics::SubSystem;

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
  ComPWA::Logging log("", "debug");

  // Construct HelicityKinematics from XML tree
  boost::property_tree::ptree tr;
  std::stringstream modelStream;
  // Construct particle list from XML tree
  modelStream << HelicityTestParticles;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, tr);

  // Construct HelicityKinematics by hand
  std::vector<int> finalState, initialState;
  initialState.push_back(421);
  finalState.push_back(310);
  finalState.push_back(-321);
  finalState.push_back(321);
  auto kin =
      std::make_shared<HelicityKinematics>(partL, initialState, finalState);

  // Generate phsp sample
  std::shared_ptr<ComPWA::Generator> gen(new ComPWA::Tools::RootGenerator(
      kin->getParticleStateTransitionKinematicsInfo(), 123));

  std::shared_ptr<ComPWA::Data::DataSet> sample(
      ComPWA::Tools::generatePhsp(20, gen));

  bool useDerivedMassSq = false;

  Event ev;
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

  double m1 = FindParticle(partL, finalState.at(0)).GetMass();
  double m2 = FindParticle(partL, finalState.at(1)).GetMass();
  double m3 = FindParticle(partL, finalState.at(2)).GetMass();
  double sqrtS = FindParticle(partL, initialState.at(0)).GetMass();

  // Define SubSystems that we want to test. The systems denoted my *_CP are
  // the CP conjugate systems the are used in the relation between D0 amplitude
  // A and D0bar amplitude Abar:
  // A(m_12^2,m_13^2) = Abar(m_13^2, m_12^2) -> A(sys12) = Aber(sys12_CP)
  // This is very specific to this decay.

  unsigned int pos_sys12(kin->addSubSystem({0}, {1}, {2}, {}));
  SubSystem sys12(kin->subSystem(pos_sys12));
  unsigned int pos_sys12_CP(kin->addSubSystem({0}, {2}, {1}, {}));
  SubSystem sys12_CP(kin->subSystem(pos_sys12_CP));

  unsigned int pos_sys13(kin->addSubSystem({2}, {0}, {1}, {}));
  SubSystem sys13(kin->subSystem(pos_sys13));
  unsigned int pos_sys13_CP(kin->addSubSystem({1}, {0}, {2}, {}));
  SubSystem sys13_CP(kin->subSystem(pos_sys13_CP));

  unsigned int pos_sys23(kin->addSubSystem({2}, {1}, {0}, {}));
  SubSystem sys23(kin->subSystem(pos_sys23));
  unsigned int pos_sys23_CP(kin->addSubSystem({1}, {2}, {0}, {}));
  SubSystem sys23_CP(kin->subSystem(pos_sys23_CP));

  LOG(INFO) << "Loop over phsp events....";
  for (auto i : sample->getEventList()) {
    // Calculate masses from FourMomentum to make sure that the correct masses
    // are used for the calculation of the helicity angle
    //    BOOST_CHECK_EQUAL((float)m1,
    //                      (float)i.getParticle(0).GetFourMomentum().GetInvMass());
    //    BOOST_CHECK_EQUAL((float)m2,
    //                      (float)i.getParticle(1).GetFourMomentum().GetInvMass());
    //    BOOST_CHECK_EQUAL((float)m3,
    //                      (float)i.getParticle(2).GetFourMomentum().GetInvMass());
    //    BOOST_CHECK_EQUAL(sqrtS, (i.getParticle(0).GetFourMomentum() +
    //                              i.getParticle(1).GetFourMomentum() +
    //                              i.getParticle(2).GetFourMomentum())
    //                                 .GetInvMass());

    double m23sq =
        (i.ParticleList[1].fourMomentum() + i.ParticleList[2].fourMomentum())
            .invMassSq();
    double m13sq =
        (i.ParticleList[0].fourMomentum() + i.ParticleList[2].fourMomentum())
            .invMassSq();
    double m12sq;
    if (useDerivedMassSq)
      m12sq = (sqrtS * sqrtS + m1 * m1 + m2 * m2 + m3 * m3 - m23sq - m13sq);
    else
      m12sq =
          (i.ParticleList[0].fourMomentum() + i.ParticleList[1].fourMomentum())
              .invMassSq();

    double RelativeTolerance(1e-6);
    //------------ Restframe (12) -------------
    // Angle in the rest frame of (12) between (1) and (3)
    double cosTheta12_13 = kin->helicityAngle(sqrtS, m1, m2, m3, m12sq, m13sq);
    double cosTheta12_13_2 =
        kin->helicityAngle(sqrtS, m2, m1, m3, m12sq, m23sq);
    // Equal to the same angle calculated from m23sq
    BOOST_CHECK_CLOSE(cosTheta12_13, -1 * cosTheta12_13_2, RelativeTolerance);

    // Angle in the rest frame of (12) between (2) and (3)
    double cosTheta12_23 = kin->helicityAngle(sqrtS, m2, m1, m3, m12sq, m23sq);
    double cosTheta12_23_2 =
        kin->helicityAngle(sqrtS, m1, m2, m3, m12sq, m13sq);
    // Equal to the same angle calculated from m13sq
    BOOST_CHECK_CLOSE(cosTheta12_23, -1 * cosTheta12_23_2, RelativeTolerance);

    BOOST_CHECK_CLOSE(cosTheta12_13, -1 * cosTheta12_23, RelativeTolerance);
    BOOST_CHECK_CLOSE(cosTheta12_13, -1 * cosTheta12_23, RelativeTolerance);
    BOOST_CHECK_CLOSE(cosTheta12_13_2, -1 * cosTheta12_23_2, RelativeTolerance);

    //------------ Restframe (13) -------------
    // Angle in the rest frame of (13) between (1) and (2)
    double cosTheta13_12 = kin->helicityAngle(sqrtS, m1, m3, m2, m13sq, m12sq);
    double cosTheta13_12_2 =
        kin->helicityAngle(sqrtS, m3, m1, m2, m13sq, m23sq);
    BOOST_CHECK_CLOSE(cosTheta13_12, -1 * cosTheta13_12_2, RelativeTolerance);
    // Angle in the rest frame of (13) between (3) and (2)
    double cosTheta13_23 = kin->helicityAngle(sqrtS, m3, m1, m2, m13sq, m23sq);
    double cosTheta13_23_2 =
        kin->helicityAngle(sqrtS, m1, m3, m2, m13sq, m12sq);
    BOOST_CHECK_CLOSE(cosTheta13_23, -1 * cosTheta13_23_2, RelativeTolerance);
    BOOST_CHECK_CLOSE(cosTheta13_12, -1 * cosTheta13_23, RelativeTolerance);
    BOOST_CHECK_CLOSE(cosTheta13_12_2, -1 * cosTheta13_23_2, RelativeTolerance);

    //------------ Restframe (23) -------------
    // Angle in the rest frame of (23) between (2) and (1)
    double cosTheta23_12 = kin->helicityAngle(sqrtS, m2, m3, m1, m23sq, m12sq);
    double cosTheta23_12_2 =
        kin->helicityAngle(sqrtS, m3, m2, m1, m23sq, m13sq);
    BOOST_CHECK_CLOSE(cosTheta23_12, -1 * cosTheta23_12_2, RelativeTolerance);

    // Angle in the rest frame of (23) between (3) and (1)
    double cosTheta23_13 = kin->helicityAngle(sqrtS, m3, m2, m1, m23sq, m13sq);
    double cosTheta23_13_2 =
        kin->helicityAngle(sqrtS, m2, m3, m1, m23sq, m12sq);
    BOOST_CHECK_CLOSE(cosTheta23_13, -1 * cosTheta23_13_2, RelativeTolerance);

    BOOST_CHECK_CLOSE(cosTheta23_12, -1 * cosTheta23_13, RelativeTolerance);

    //------------- Test of HelicityKinematics -----------------
    // Check if calculation of helicity angles corresponds to the previously
    // calculated values

    DataPoint p12, p12_CP, p13, p13_CP, p23, p23_CP;

    LOG(DEBUG) << "-------- NEW EVENT ----------";
    kin->convert(i, p12, sys12);
    kin->convert(i, p12_CP, sys12_CP);
    // BOOST_CHECK_EQUAL((float)p12.value(1), (float)cosTheta12_23);
    BOOST_CHECK(ComPWA::equal(std::cos(p12.KinematicVariableList[1]),
                              cosTheta12_23, 1000));

    LOG(DEBUG) << "-------- (12) ----------";
    LOG(DEBUG) << sys12 << " : " << p12;
    LOG(DEBUG) << sys12_CP << " : " << p12_CP;
    LOG(DEBUG) << "cosTheta12 "
               << kin->helicityAngle(sqrtS, m2, m1, m3, m12sq, m23sq) << " CP: "
               << kin->helicityAngle(sqrtS, m3, m1, m2, m13sq, m23sq);

    kin->convert(i, p13, sys13);
    kin->convert(i, p13_CP, sys13_CP);
    BOOST_CHECK_CLOSE(std::cos(p13.KinematicVariableList[1]), cosTheta13_12,
                      RelativeTolerance);

    BOOST_CHECK_CLOSE(std::cos(p13.KinematicVariableList[1]),
                      -1 * std::cos(p12_CP.KinematicVariableList[1]),
                      RelativeTolerance);
    BOOST_CHECK_CLOSE(std::cos(p12.KinematicVariableList[1]),
                      -1 * std::cos(p13_CP.KinematicVariableList[1]),
                      RelativeTolerance);

    LOG(DEBUG) << "-------- (13) ----------";
    LOG(DEBUG) << sys13 << " : " << p13;
    LOG(DEBUG) << sys13_CP << " : " << p13_CP;
    LOG(DEBUG) << "cosTheta13 "
               << kin->helicityAngle(sqrtS, m1, m3, m2, m13sq, m12sq) << " CP: "
               << kin->helicityAngle(sqrtS, m1, m2, m3, m12sq, m13sq);

    kin->convert(i, p23, sys23);
    kin->convert(i, p23_CP, sys23_CP);
    BOOST_CHECK_CLOSE(std::cos(p23.KinematicVariableList[1]), cosTheta23_12,
                      RelativeTolerance);
    BOOST_CHECK_CLOSE(std::cos(p23.KinematicVariableList[1]),
                      -1 * std::cos(p23_CP.KinematicVariableList[1]),
                      RelativeTolerance);

    LOG(DEBUG) << "-------- (23) ----------";
    LOG(DEBUG) << sys23 << " : " << p23;
    LOG(DEBUG) << sys23_CP << " : " << p23_CP;
    LOG(DEBUG) << "cosTheta23 "
               << kin->helicityAngle(sqrtS, m2, m3, m1, m23sq, m12sq) << " CP: "
               << kin->helicityAngle(sqrtS, m3, m2, m1, m23sq, m13sq);
  }
}

BOOST_AUTO_TEST_SUITE_END()
