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

#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/qft++/Vector4.h"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

HelicityKinematics::HelicityKinematics(std::shared_ptr<PartList> partL,
                                       std::vector<pid> initialState,
                                       std::vector<pid> finalState)
    : Kinematics(initialState, finalState), _partList(partL) {

  // Currently we can only handle a particle decay.
  // Note: Whats different for particle scattering?
  assert(_initialState.size() == 1);

  // Create title
  std::stringstream stream;
  stream << "( ";
  for (auto i : _initialState)
    stream << FindParticle(partL, i).GetName() << " ";
  stream << ")->( ";
  for (auto i : _finalState)
    stream << FindParticle(partL, i).GetName() << " ";
  stream << ")";

  LOG(info) << "HelicityKinematics::HelicityKinematics() | Initialize reaction "
            << stream.str();
  return;
}

HelicityKinematics::HelicityKinematics(std::shared_ptr<PartList> partL,
                                       boost::property_tree::ptree pt)
    : _partList(partL) {

  auto initialS = pt.get_child("InitialState");
  _initialState = std::vector<int>(initialS.size());
  for (auto i : initialS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = partL->find(name)->second;
    int pos = i.second.get<int>("<xmlattr>.Id");
    _initialState.at(pos) = partP.GetId();
  }
  // Currently we can only handle a particle decay.
  // Note: Whats different for particle scattering?
  assert(_initialState.size() == 1);

  auto finalS = pt.get_child("FinalState");
  _finalState = std::vector<int>(finalS.size());
  for (auto i : finalS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = partL->find(name)->second;
    int pos = i.second.get<int>("<xmlattr>.Id");
    _finalState.at(pos) = partP.GetId();
  }

  auto phspVal = pt.get_optional<double>("PhspVolume");
  if (phspVal) {
    SetPhspVolume(phspVal.get());
  }

  // Creating unique title
  std::stringstream stream;
  stream << "( ";
  for (auto i : _initialState)
    stream << FindParticle(partL, i).GetName() << " ";
  stream << ")->( ";
  for (auto i : _finalState)
    stream << FindParticle(partL, i).GetName() << " ";
  stream << ")";

  LOG(info) << "HelicityKinematics::HelicityKinematics() | Initialize reaction "
            << stream.str();
  return;
}

bool HelicityKinematics::IsWithinPhsp(const dataPoint &point) const {
  int subSystemID = 0;
  int pos = 0;
  while ((pos + 2) < point.Size()) {
    auto invMassBounds = GetInvMassBounds(GetSubSystem(subSystemID));
    if (point.GetValue(pos) < invMassBounds.first ||
        point.GetValue(pos) > invMassBounds.second)
      return false;
    if (point.GetValue(pos + 1) < -1 || point.GetValue(pos + 1) > +1)
      return false;
    if (point.GetValue(pos + 2) < 0 || point.GetValue(pos + 2) > 2 * M_PI)
      return false;

    pos += 3;
    subSystemID++;
  }

  return true;
}

void HelicityKinematics::EventToDataPoint(const Event &event,
                                          dataPoint &point) const {
  assert(_listSubSystem.size() == _invMassBounds.size());

  for (int i = 0; i < _listSubSystem.size(); i++)
    EventToDataPoint(event, point, _listSubSystem.at(i), _invMassBounds.at(i));
  return;
}

void HelicityKinematics::EventToDataPoint(const Event &event, dataPoint &point,
                                          const SubSystem &sys) const {
  auto massLimits = GetInvMassBounds(sys);
  EventToDataPoint(event, point, sys, massLimits);
}

double HelicityKinematics::HelicityAngle(double M, double m, double m2,
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

void HelicityKinematics::EventToDataPoint(
    const Event &event, dataPoint &point, const SubSystem &sys,
    const std::pair<double, double> limits) const {

  assert(sys.GetFinalStates().size() == 2 &&
         "HelicityKinematics::EventToDataPoint() | More then two particles.");

  FourMomentum cms;
  for (auto s : sys.GetRecoilState())
    cms += event.GetParticle(s).GetFourMomentum();

  FourMomentum finalA, finalB;
  for (auto s : sys.GetFinalStates().at(0))
    finalA += event.GetParticle(s).GetFourMomentum();

  for (auto s : sys.GetFinalStates().at(1))
    finalB += event.GetParticle(s).GetFourMomentum();

  // Four momentum of the decaying resonance
  FourMomentum resP4 = finalA + finalB;
  double mSq = resP4.GetInvMassSq();
  
  // Calculate sum of final states four momenta
  cms += resP4;

  if (mSq <= limits.first) {
    // We allow for a deviation from the limits of 10 times the numerical
    // precision
    if (ComPWA::equal(mSq, limits.first, 10))
      mSq = limits.first;
    else
      throw BeyondPhsp("HelicityKinematics::EventToDataPoint() |"
                       " Point beypond phase space boundaries!");
  }
  if (mSq >= limits.second) {
    // We allow for a deviation from the limits of 10 times the numerical
    // precision
    if (ComPWA::equal(mSq, limits.second, 10))
      mSq = limits.second;
    else
      throw BeyondPhsp("HelicityKinematics::EventToDataPoint() |"
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
    throw BeyondPhsp("HelicityKinematics::EventToDataPoint() |"
                     " Point beypond phase space boundaries!");
  }

  point.GetPoint().push_back(mSq);
  point.GetPoint().push_back(cosTheta);
  point.GetPoint().push_back(phi);
}

const std::pair<double, double> &
HelicityKinematics::GetInvMassBounds(const SubSystem &sys) const {
  return GetInvMassBounds(
      const_cast<HelicityKinematics *>(this)->GetDataID(sys));
}

const std::pair<double, double> &
HelicityKinematics::GetInvMassBounds(int sysID) const {
  return _invMassBounds.at(sysID);
}

std::pair<double, double>
HelicityKinematics::CalculateInvMassBounds(const SubSystem &sys) const {

  /// We use the formulae from (PDG2016 Kinematics Fig.47.3). I hope the
  /// generalization to n-body decays is correct.
  std::pair<double, double> lim(
      0, FindParticle(_partList, _initialState.at(0)).GetMass());
  // Sum up masses of all final state particles
  for (auto j : sys.GetFinalStates())
    for (auto i : j)
      lim.first += FindParticle(_partList, _finalState.at(i)).GetMass();
  lim.first *= lim.first;

  for (auto i : sys.GetRecoilState())
    lim.second -= FindParticle(_partList, _finalState.at(i)).GetMass();
  lim.second *= lim.second;

  return lim;
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
