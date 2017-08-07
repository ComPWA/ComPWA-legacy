// Define Boost test module
#define BOOST_TEST_MODULE HelicityFormalism

#include <vector>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "Core/PhysConst.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Logging.hpp"
#include "Core/Particle.hpp"
#include "Tools/RunManager.hpp"
#include "Tools/RootGenerator.hpp"
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

using namespace ComPWA;
using namespace ComPWA::Physics::HelicityFormalism;

// Define Boost test suite (no idea what's the difference to TEST_MODULE)
BOOST_AUTO_TEST_SUITE(HelicityFormalism)

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
  ComPWA::Logging log("", boost::log::trivial::severity_level::debug);

  // Construct HelicityKinematics from XML tree
  boost::property_tree::ptree tr;
  boost::property_tree::xml_parser::read_xml("AmpModel-input.xml", tr);

  ComPWA::PhysConst::CreateInstance(tr);

  // Construct HelicityKinematics by hand
  std::vector<int> finalState, initialState;
  initialState.push_back(421);
  finalState.push_back(310);
  finalState.push_back(-321);
  finalState.push_back(321);
  auto kin = std::make_shared<HelicityKinematics>(initialState, finalState);

  // Generate phsp sample
  std::shared_ptr<ComPWA::Generator> gen(new ComPWA::Tools::RootGenerator(
      kin->GetInitialState(),
      kin->GetFinalState(), 123));
  std::shared_ptr<ComPWA::DataReader::Data> sample(
      new ComPWA::DataReader::RootReader());

  ComPWA::RunManager r;
  r.SetGenerator(gen);
  r.SetPhspSample(sample);
  r.GeneratePhsp(20);

  bool useDerivedMassSq = false;

  Event ev;
  /* We add an event of D0->KsK-K+ from data. Since the decay KS->pipi is not
   * constraint to the KS mass we use the relation:
   * (sqrtS * sqrtS + m1 * m1 + m2 * m2 + m3 * m3 - m23sq - m13sq)
   * to calculated the third invariant mass. The results for the sub systems
   * listed below are:
   * m23sq=1.4014 m13sq=1.52861 m12sq=1.28267
   * cosTheta12_23=0.151776 cosTheta12_CP=-0.178456
   * cosTheta13_12=0.178456 cosTheta13_CP=-0.151776
   * cosTheta23_12=0.318779 cosTheta23_CP=-0.318779
   * Angles are symmetric! Set @useDerivedMassSq = true.

   * In case all invariant masses are calculated independently
   * (@useDerivedMassSq = false)the angles are not symmetric anymore.
   * In this case a couple of tests are supposed to fail.
   * m23sq=1.4014 m13sq=1.52861 m12sq=1.26798
   * cosTheta12_23=0.171554 cosTheta12_CP=-0.178456
   * cosTheta13_12=0.219724 cosTheta13_CP=-0.13337
   * cosTheta23_12=0.356843 cosTheta23_CP=-0.318779
   */
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

  double m1 = PhysConst::Instance()->FindParticle(finalState.at(0)).GetMass();
  double m2 = PhysConst::Instance()->FindParticle(finalState.at(1)).GetMass();
  double m3 = PhysConst::Instance()->FindParticle(finalState.at(2)).GetMass();
  double sqrtS =
      PhysConst::Instance()->FindParticle(initialState.at(0)).GetMass();

  // Define SubSystems that we want to test. The systems denoted my *_CP are
  // the CP conjugate systems the are used in the relation between D0 amplitude
  // A and D0bar amplitude Abar:
  // A(m_12^2,m_13^2) = Abar(m_13^2, m_12^2) -> A(sys12) = Aber(sys12_CP)
  // This is very specific to this decay.
  SubSystem sys12(std::vector<int>{2}, std::vector<int>{1},
                  std::vector<int>{0});
  SubSystem sys12_CP(std::vector<int>{1}, std::vector<int>{2},
                     std::vector<int>{0});

  SubSystem sys13(std::vector<int>{1}, std::vector<int>{0},
                  std::vector<int>{2});
  SubSystem sys13_CP(std::vector<int>{2}, std::vector<int>{0},
                     std::vector<int>{1});

  SubSystem sys23(std::vector<int>{0}, std::vector<int>{1},
                  std::vector<int>{2});
  SubSystem sys23_CP(std::vector<int>{0}, std::vector<int>{2},
                     std::vector<int>{1});

  LOG(info) << "Loop over phsp events....";
  for (auto i : sample->GetEvents()) {
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

    double m23sq = (i.GetParticle(1).GetFourMomentum() +
                    i.GetParticle(2).GetFourMomentum())
                       .GetInvMassSq();
    double m13sq = (i.GetParticle(0).GetFourMomentum() +
                    i.GetParticle(2).GetFourMomentum())
                       .GetInvMassSq();
    double m12sq;
    if (useDerivedMassSq)
      m12sq = (sqrtS * sqrtS + m1 * m1 + m2 * m2 + m3 * m3 - m23sq - m13sq);
    else
      m12sq = (i.GetParticle(0).GetFourMomentum() +
               i.GetParticle(1).GetFourMomentum())
                  .GetInvMassSq();

    //------------ Restframe (12) -------------
    // Angle in the rest frame of (12) between (1) and (3)
    double cosTheta12_13 = kin->HelicityAngle(sqrtS, m1, m2, m3, m12sq, m13sq);
    double cosTheta12_13_2 =
        kin->HelicityAngle(sqrtS, m2, m1, m3, m12sq, m23sq);
    // Equal to the same angle calculated from m23sq
    BOOST_CHECK_EQUAL((float)cosTheta12_13, (-1) * (float)(cosTheta12_13_2));

    // Angle in the rest frame of (12) between (2) and (3)
    double cosTheta12_23 = kin->HelicityAngle(sqrtS, m2, m1, m3, m12sq, m23sq);
    double cosTheta12_23_2 =
        kin->HelicityAngle(sqrtS, m1, m2, m3, m12sq, m13sq);
    // Equal to the same angle calculated from m13sq
    BOOST_CHECK_EQUAL((float)cosTheta12_23, (-1) * (float)(cosTheta12_23_2));

    BOOST_CHECK_EQUAL((float)cosTheta12_13, (-1) * (float)(cosTheta12_23));
    BOOST_CHECK_EQUAL((float)cosTheta12_13, (-1) * ((float)cosTheta12_23));
    BOOST_CHECK_EQUAL((float)cosTheta12_13_2, (-1) * ((float)cosTheta12_23_2));

    //------------ Restframe (13) -------------
    // Angle in the rest frame of (13) between (1) and (2)
    double cosTheta13_12 = kin->HelicityAngle(sqrtS, m1, m3, m2, m13sq, m12sq);
    double cosTheta13_12_2 =
        kin->HelicityAngle(sqrtS, m3, m1, m2, m13sq, m23sq);
    BOOST_CHECK_EQUAL((float)cosTheta13_12, (-1) * (float)(cosTheta13_12_2));
    // Angle in the rest frame of (13) between (3) and (2)
    double cosTheta13_23 = kin->HelicityAngle(sqrtS, m3, m1, m2, m13sq, m23sq);
    double cosTheta13_23_2 =
        kin->HelicityAngle(sqrtS, m1, m3, m2, m13sq, m12sq);
    BOOST_CHECK_EQUAL((float)cosTheta13_23, (-1) * ((float)cosTheta13_23_2));
    BOOST_CHECK_EQUAL((float)cosTheta13_12, (-1) * ((float)cosTheta13_23));
    BOOST_CHECK_EQUAL((float)cosTheta13_12_2, (-1) * ((float)cosTheta13_23_2));

    //------------ Restframe (23) -------------
    // Angle in the rest frame of (23) between (2) and (1)
    double cosTheta23_12 = kin->HelicityAngle(sqrtS, m2, m3, m1, m23sq, m12sq);
    double cosTheta23_12_2 =
        kin->HelicityAngle(sqrtS, m3, m2, m1, m23sq, m13sq);
    BOOST_CHECK_EQUAL((float)cosTheta23_12, (float)((-1) * cosTheta23_12_2));

    // Angle in the rest frame of (23) between (3) and (1)
    double cosTheta23_13 = kin->HelicityAngle(sqrtS, m3, m2, m1, m23sq, m13sq);
    double cosTheta23_13_2 =
        kin->HelicityAngle(sqrtS, m2, m3, m1, m23sq, m12sq);
    BOOST_CHECK_EQUAL((float)cosTheta23_13, (-1) * ((float)cosTheta23_13_2));

    BOOST_CHECK_EQUAL((float)cosTheta23_12, (-1) * ((float)cosTheta23_13));

    //------------- Test of HelicityKinematics -----------------
    // Check if calculation of helicity angles correspongs the the previously
    // calculated values

    dataPoint p12, p12_CP, p13, p13_CP, p23, p23_CP;

    LOG(debug) << "-------- NEW EVENT ----------";
    kin->EventToDataPoint(i, p12, sys12);
    kin->EventToDataPoint(i, p12_CP, sys12_CP);
    BOOST_CHECK_EQUAL((float)p12.GetValue(1), (float)cosTheta12_23);

    LOG(debug) << "-------- (12) ----------";
    LOG(debug) << sys12.to_string() << " : " << p12;
    LOG(debug) << sys12_CP.to_string() << " : " << p12_CP;
    LOG(debug) << "cosTheta12 "
               << kin->HelicityAngle(sqrtS, m2, m1, m3, m12sq, m23sq) << " CP: "
               << kin->HelicityAngle(sqrtS, m3, m1, m2, m13sq, m23sq);

    kin->EventToDataPoint(i, p13, sys13);
    kin->EventToDataPoint(i, p13_CP, sys13_CP);
    BOOST_CHECK_EQUAL((float)p13.GetValue(1), (float)cosTheta13_12);

    BOOST_CHECK_EQUAL((float)p13.GetValue(1), (-1) * (float)p12_CP.GetValue(1));
    BOOST_CHECK_EQUAL((float)p12.GetValue(1), (-1) * (float)p13_CP.GetValue(1));

    LOG(debug) << "-------- (13) ----------";
    LOG(debug) << sys13.to_string() << " : " << p13;
    LOG(debug) << sys13_CP.to_string() << " : " << p13_CP;
    LOG(debug) << "cosTheta13 "
               << kin->HelicityAngle(sqrtS, m1, m3, m2, m13sq, m12sq) << " CP: "
               << kin->HelicityAngle(sqrtS, m1, m2, m3, m12sq, m13sq);

    kin->EventToDataPoint(i, p23, sys23);
    kin->EventToDataPoint(i, p23_CP, sys23_CP);
    BOOST_CHECK_EQUAL((float)p23.GetValue(1), (float)cosTheta23_12);

    BOOST_CHECK_EQUAL((float)p23.GetValue(1),
                      (-1) * ((float)p23_CP.GetValue(1)));

    LOG(debug) << "-------- (23) ----------";
    LOG(debug) << sys23.to_string() << " : " << p23;
    LOG(debug) << sys23_CP.to_string() << " : " << p23_CP;
    LOG(debug) << "cosTheta23 "
               << kin->HelicityAngle(sqrtS, m2, m3, m1, m23sq, m12sq) << " CP: "
               << kin->HelicityAngle(sqrtS, m3, m2, m1, m23sq, m13sq);
  }
}

BOOST_AUTO_TEST_SUITE_END()
