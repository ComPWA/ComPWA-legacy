// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <numeric>
#include <cmath>
#include <algorithm>

#include "Core/Event.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Particle.hpp"
#include "Core/Properties.hpp"

#include "Physics/EvtGen/DalitzKinematics.hpp"
#include "ThirdParty/qft++/include/qft++/Vector4.h"

namespace ComPWA {
namespace Physics {
namespace EvtGenIF {

DalitzKinematics::DalitzKinematics(std::shared_ptr<PartList> partL,
                                       std::vector<pid> initialState,
                                       std::vector<pid> finalState,
                                       ComPWA::FourMomentum cmsP4)
    : Kinematics(initialState, finalState, cmsP4), ParticleList(partL),
	  _M(FindParticle(partL, InitialState.at(0)).GetMass()) {

  // If the cms four-momentum is not set of set it here
  if (InitialStateP4 == FourMomentum(0, 0, 0, 0) && InitialState.size() == 1) {
      double sqrtS = FindParticle(partL, InitialState.at(0)).GetMass();
      InitialStateP4 = ComPWA::FourMomentum(0, 0, 0, sqrtS);
  }

  // Make sure cms momentum is set
  if (InitialStateP4 == FourMomentum(0, 0, 0, 0))
    assert(false);

  // Create title
  std::stringstream stream;
  stream << "( ";
  for (auto i : InitialState)
    stream << FindParticle(partL, i).name() << " ";
  stream << ")->( ";
  for (auto i : InitialState)
    stream << FindParticle(partL, i).name() << " ";
  stream << ")";

  LOG(INFO) << "DalitzKinematics::DalitzKinematics() | Initialize reaction "
            << stream.str();
  return;
}

DalitzKinematics::DalitzKinematics(std::shared_ptr<PartList> partL,
                                       boost::property_tree::ptree pt)
    : ParticleList(partL),_M(FindParticle(partL, InitialState.at(0)).GetMass()) {

  auto initialS = pt.get_child("InitialState");
  InitialState = std::vector<int>(initialS.size());
  for (auto i : initialS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = partL->find(name)->second;
    int pos = i.second.get<int>("<xmlattr>.Id");
    InitialState.at(pos) = partP.GetId();
  }
  if (InitialState.size() == 1) {
    double sqrtS = FindParticle(partL, InitialState.at(0)).GetMass();
    InitialStateP4 = ComPWA::FourMomentum(0, 0, 0, sqrtS);
  } else {
    // If the initial state has more than one particle we require that the
    // initial four momentum is specified
    InitialStateP4 = FourMomentumFactory(pt.get_child("InitialFourMomentum"));
  }

  // Make sure cms momentum is set
  if (InitialStateP4 == FourMomentum(0, 0, 0, 0))
    assert(false);

  auto finalS = pt.get_child("FinalState");
  FinalState = std::vector<int>(finalS.size());
  for (auto i : finalS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = partL->find(name)->second;
    int pos = i.second.get<int>("<xmlattr>.Id");
    FinalState.at(pos) = partP.GetId();
  }

  auto phspVal = pt.get_optional<double>("PhspVolume");
  if (phspVal) {
    setPhspVolume(phspVal.get());
  }

  // Creating unique title
  std::stringstream stream;
  stream << "( ";
  for (auto i : InitialState)
    stream << FindParticle(partL, i).name() << " ";
  stream << ")->( ";
  for (auto i : FinalState)
    stream << FindParticle(partL, i).name() << " ";
  stream << ")";

  LOG(INFO) << "DalitzKinematics::DalitzKinematics() | Initialize reaction "
            << stream.str();
  return;
}

bool DalitzKinematics::isWithinPhsp(const DataPoint &point) const {
  int subSystemID = 0;
  int pos = 0;

  double mA = point.value(0);
  double mB = point.value(1);
  double mC = point.value(2);
  double qAB = point.value(3);
  double qBC = point.value(4);
  double qCA = point.value(5);

  double s2 = (qAB+qBC+qCA - mA*mA - mB*mB - mC*mC);

  if(s2<(_M*_M)) return true;

  return false;
}

void DalitzKinematics::convert(const Event &event,
                                          DataPoint &point) const {

	double mA,mB,mC,qAB,qBC,qCA;

	  mA = event.particle(0).mass();
	  mB = event.particle(1).mass();
	  mC = event.particle(2).mass();
	  qAB = FourMomentum::invariantMass(event.particle(0).fourMomentum(), event.particle(1).fourMomentum());
	  qBC = FourMomentum::invariantMass(event.particle(1).fourMomentum(), event.particle(2).fourMomentum());
	  qCA = FourMomentum::invariantMass(event.particle(2).fourMomentum(), event.particle(0).fourMomentum());

	  /*FourMomentum cms;
	  for (auto s : sys.GetRecoilState())
	    cms += event.particle(s).fourMomentum();

	  FourMomentum finalA, finalB;
	  for (auto s : sys.GetFinalStates().at(0))
	    finalA += event.particle(s).fourMomentum();

  assert(sys.GetFinalStates().size() == 2 &&
         "DalitzKinematics::convert() | More then two particles.");


  for (auto s : sys.GetFinalStates().at(1))
    finalB += event.particle(s).fourMomentum();

  // Four momentum of the decaying resonance
  FourMomentum resP4 = finalA + finalB;
  double mSq = resP4.invMassSq();

  // Calculate sum of final states four momenta
  cms += resP4;

  if (mSq <= limits.first) {
    // We allow for a deviation from the limits of 10 times the numerical
    // precision
    if (ComPWA::equal(mSq, limits.first, 10))
      mSq = limits.first;
    else
      throw BeyondPhsp("HelicityKinematics::convert() |"
                       " Point beypond phase space boundaries!");
  }
  if (mSq >= limits.second) {
    // We allow for a deviation from the limits of 10 times the numerical
    // precision
    if (ComPWA::equal(mSq, limits.second, 10))
      mSq = limits.second;
    else
      throw BeyondPhsp("HelicityKinematics::convert() |"
                       " Point beypond phase space boundaries!");
  }

  // When using finalB instead of finalA here, the WignerD changes sign. In
  // the end this does not matter
  QFT::Vector4<double> p4QftCms(cms);
  QFT::Vector4<double> p4QftResonance(resP4);
  QFT::Vector4<double> p4QftDaughter(finalB);

  // Boost the four momentum of the decaying resonance to total CMS
  p4QftResonance.Boost(p4QftCms);
  // Boost the four momentum of one daughter particle to CMS of the resonance
  p4QftDaughter.Boost(p4QftResonance);

  // Calculate the angles between recoil system and final state.
  // Use an Euler rotation of the coordinate system (wrong?)
  //   p4QftDaughter.Rotate(p4QftResonance.Phi(), p4QftResonance.Theta(),
  //                       (-1) * p4QftResonance.Phi());
  p4QftDaughter.RotateZ((-1) * p4QftResonance.Phi());
  p4QftDaughter.RotateY((-1) * p4QftResonance.Theta());

  double cosTheta = p4QftDaughter.CosTheta();
  double phi = p4QftDaughter.Phi();

  //  double cc;
  //  if (sys.GetRecoilState().size() == 1 &&
  //      sys.GetFinalStates().at(0).size() == 1 &&
  //      sys.GetFinalStates().at(1).size() == 1) {
  //    double invMassSqA = mSq;
  //    double invMassSqB = (recoilP4 + finalA).GetInvMassSq();
  //    auto mspec = PhysConst::Instance()
  //                     ->FindParticle(_finalState.at(sys.GetRecoilState().at(0)))
  //                     .GetMass();
  //    auto ma =
  //        PhysConst::Instance()
  //            ->FindParticle(_finalState.at(sys.GetFinalStates().at(0).at(0)))
  //            .GetMass();
  //    auto mb =
  //        PhysConst::Instance()
  //            ->FindParticle(_finalState.at(sys.GetFinalStates().at(1).at(0)))
  //            .GetMass();
  //    auto M =
  //    PhysConst::Instance()->FindParticle(_initialState.at(0)).GetMass();
  //
  //    cc = HelicityAngle(M, ma, mb, mspec, invMassSqA, invMassSqB);
  //    std::cout << sys << std::endl;
  //    std::cout << _initialState.at(0) << "/ (" <<
  // sys.GetFinalStates().at(0).at(0)
  //              << sys.GetFinalStates().at(1).at(0) << ") angle ("<<
  // sys.GetFinalStates().at(0).at(0)
  //              << sys.GetRecoilState().at(0) << ") - "
  //              << " (ma=" << ma<<" mb="<<mb<<" mSpec="<<mspec<<"
  // mABSq="<<invMassSqA << " mASpecSq=" << invMassSqB << ") = " << cc
  //              << " " << cosTheta << std::endl;
  //
  //  } else {
  //    cc = 1.0;
  //    phi = 0.0;
  //  }
  //  cosTheta = cc;

  //   Check if values are within allowed range.
  if (cosTheta > 1 || cosTheta < -1 || phi > M_PI || phi < (-1) * M_PI ||
      std::isnan(cosTheta) || std::isnan(phi)) {
    throw BeyondPhsp("HelicityKinematics::convert() |"
                     " Point beypond phase space boundaries!");
  }*/

  //point.values().push_back(mSq);
  //point.values().push_back(cosTheta);
  //point.values().push_back(phi);
  point.values().push_back(mA);
  point.values().push_back(mB);
  point.values().push_back(mC);
  point.values().push_back(qAB);
  point.values().push_back(qBC);
  point.values().push_back(qCA);


  return;
}

double DalitzKinematics::helicityAngle(double M, double m, double m2,
                                         double mSpec, double invMassSqA,
                                         double invMassSqB) const {
  // Calculate energy and momentum of m1/m2 in the invMassSqA rest frame
  double eCms = (invMassSqA + m * m - m2 * m2) / (2 * sqrt(invMassSqA));
  double pCms = eCms * eCms - m * m;
  // Calculate energy and momentum of mSpec in the invMassSqA rest frame
  double eSpecCms =
      (M * M - mSpec * mSpec - invMassSqA) / (2 * sqrt(invMassSqA));
  double pSpecCms = eSpecCms * eSpecCms - mSpec * mSpec;
  double cosAngle =
      -(invMassSqB - m * m - mSpec * mSpec - 2 * eCms * eSpecCms) /
      (2 * sqrt(pCms * pSpecCms));

  //  if( cosAngle>1 || cosAngle<-1 ){
  //      throw BeyondPhsp("DalitzKinematics::helicityAngle() | "
  //              "scattering angle out of range! Datapoint beyond phsp? angle="
  //              +std::to_string((long double) cosAngle)
  //      +" M="+std::to_string((long double) M)
  //      +" m="+std::to_string((long double) m)
  //      +" m2="+std::to_string((long double) m2)
  //      +" mSpec="+std::to_string((long double) mSpec)
  //      +" mSqA="+std::to_string((long double) invMassSqA)
  //      +" mSqB="+std::to_string((long double) invMassSqB) );
  //  }
  if (cosAngle > 1 || cosAngle < -1) { // faster
    throw BeyondPhsp("DalitzKinematics::helicityAngle() | "
                     "scattering angle out of range! Datapoint beyond"
                     "phsp?");
  }
  return cosAngle;
}


int DalitzKinematics::createIndex() {

	if(VariableNames.size()) return 1;

    VariableNames.push_back("mA");
    VariableNames.push_back("mB");
    VariableNames.push_back("mC");
    VariableNames.push_back("qAB");
    VariableNames.push_back("qBC");
    VariableNames.push_back("qCA");
    VariableTitles.push_back("mA");
    VariableTitles.push_back("mB");
    VariableTitles.push_back("mC");
    VariableTitles.push_back("qAB");
    VariableTitles.push_back("qBC");
    VariableTitles.push_back("qCA");
  //}
  return 0;
}

} // ns::EvtGenIF
} // ns::Physics
} // ns::ComPWA
