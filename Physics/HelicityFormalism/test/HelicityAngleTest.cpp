#define BOOST_TEST_MODULE                                                      \
  HelicityFormalism /* this can only be define once within the same library ?! \
  */
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/PhysConst.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Logging.hpp"
#include "Core/Particle.hpp"
#include "Tools/RunManager.hpp"
#include "Tools/RootGenerator.hpp"
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

using namespace ComPWA::Physics::HelicityFormalism;

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <iostream>

using namespace ComPWA;

BOOST_AUTO_TEST_SUITE(HelicityFormalism)

//TODO: General description
BOOST_AUTO_TEST_CASE(HelicityAngleTest) {
  ComPWA::Logging log("", boost::log::trivial::severity_level::error);

  // Construct HelicityKinematics from XML tree
  boost::property_tree::ptree tr;
  boost::property_tree::xml_parser::read_xml("HelicityFormalismTest-input.xml",
                                             tr);

  ComPWA::PhysConst::CreateInstance(tr);

  // Construct HelicityKinematics by hand
  std::vector<int> finalState, initialState;
  initialState.push_back(421);
  finalState.push_back(310);
  finalState.push_back(-321);
  finalState.push_back(321);
  HelicityKinematics::CreateInstance(initialState, finalState);

  // Generate phsp sample
  std::shared_ptr<ComPWA::Generator> gen(new ComPWA::Tools::RootGenerator(123));
  std::shared_ptr<ComPWA::DataReader::Data> sample(
      new ComPWA::DataReader::RootReader());

  ComPWA::RunManager r;
  r.SetGenerator(gen);
  r.SetPhspSample(sample);
  r.GeneratePhsp(5);

  bool useDerivedMassSq = true;

  Event ev;
  /* We add an event from data. We do so since due to measurement errors
   * the symmetry between sets of helicity angles are destroyed. This is
   * especially the case if all three invariant masses are directly
   * calculated from the four-momenta and the relation:
   * (sqrtS * sqrtS + m1 * m1 + m2 * m2 + m3 * m3 - m23sq - m13sq)
   * is not used. This can bw seed from the following (hopefully) correct
   * results for the angles:
   *
   * All invariant masses calculated from four-momenta. Angles are not
   * symmetric! Set @useDerivedMassSq = false.
   * is this case a couple of tests are supposed to fail.
   * m23sq=1.4014 m13sq=1.52861 m12sq=1.26798
   * cosTheta12_23=0.171554 cosTheta12_CP=-0.178456
   * cosTheta13_12=0.219724 cosTheta13_CP=-0.13337
   * cosTheta23_12=0.356843 cosTheta23_CP=-0.318779
   *
   * The third invariant mass (m12sq) is derived from the others
   * (m23sq and m13sq). Angles are symmetric! Set @useDerivedMassSq = true.
   * m23sq=1.4014 m13sq=1.52861 m12sq=1.28267
   * cosTheta12_23=0.151776 cosTheta12_CP=-0.178456
   * cosTheta13_12=0.178456 cosTheta13_CP=-0.151776
   * cosTheta23_12=0.318779 cosTheta23_CP=-0.318779
   */
  ev.addParticle(ComPWA::Particle(
      std::array<double, 4>{{-0.00827061, -0.242581, -0.335833, 0.636104}},
      310));
  ev.addParticle(ComPWA::Particle(
      std::array<double, 4>{{-0.158637, -0.149132, 0.199913, 0.575405}}, 321));
  ev.addParticle(ComPWA::Particle(
      std::array<double, 4>{{-0.0236227, 0.453598, -0.0330521, 0.671656}},
      -321));
  sample->pushEvent(ev);

  auto kin = dynamic_cast<HelicityKinematics *>(Kinematics::Instance());

  double m1 = PhysConst::Instance()->FindParticle(finalState.at(0)).GetMass();
  double m2 = PhysConst::Instance()->FindParticle(finalState.at(1)).GetMass();
  double m3 = PhysConst::Instance()->FindParticle(finalState.at(2)).GetMass();
  double sqrtS =
      PhysConst::Instance()->FindParticle(initialState.at(0)).GetMass();

  // Define SubSystems that we want to test
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
  for (auto i : sample->getEvents()) {
    // Calculate masses from FourMomentum to make sure that the correct masses
    // are used for the calculation of the helicity angle
    BOOST_CHECK_EQUAL((float)m1,
                      (float)i.getParticle(0).GetFourMomentum().GetInvMass());
    BOOST_CHECK_EQUAL((float)m2,
                      (float)i.getParticle(1).GetFourMomentum().GetInvMass());
    BOOST_CHECK_EQUAL((float)m3,
                      (float)i.getParticle(2).GetFourMomentum().GetInvMass());
    //    double sqrtS = (i.getParticle(0).GetFourMomentum() +
    //                    i.getParticle(1).GetFourMomentum() +
    //                    i.getParticle(2).GetFourMomentum())
    //                       .GetInvMass();

    double m23sq = (i.getParticle(1).GetFourMomentum() +
                    i.getParticle(2).GetFourMomentum())
                       .GetInvMassSq();
    double m13sq = (i.getParticle(0).GetFourMomentum() +
                    i.getParticle(2).GetFourMomentum())
                       .GetInvMassSq();
    double m12sq;
    if (useDerivedMassSq)
      m12sq = (sqrtS * sqrtS + m1 * m1 + m2 * m2 + m3 * m3 - m23sq - m13sq);
    else
      m12sq = (i.getParticle(0).GetFourMomentum() +
               i.getParticle(1).GetFourMomentum())
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
