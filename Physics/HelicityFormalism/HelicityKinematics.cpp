// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <numeric>
#include <cmath>
#include <algorithm>

#include "Core/Event.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Particle.hpp"
#include "Core/PhysConst.hpp"
#include "DataReader/RootReader/RootReader.hpp"

#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/qft++/Vector4.h"
#include "Tools/RootGenerator.hpp"

#include <boost/timer.hpp>

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

HelicityKinematics::HelicityKinematics(std::vector<pid> initialState,
                                       std::vector<pid> finalState)
    : Kinematics(initialState, finalState) {
  assert(_initialState.size() == 1);
  pid idMother = _initialState.at(0);
  auto motherProp = PhysConst::Instance()->FindParticle(idMother);
  _M = motherProp.GetMass();
  _Msq = _M * _M;
  _spinM = motherProp.GetSpinQuantumNumber("Spin");

  // Creating unique title
  std::stringstream stream;
  stream << "( ";
  for (auto i : _initialState)
    stream << PhysConst::Instance()->FindParticle(i).GetName() << " ";
  stream << ")->( ";
  for (auto i : _finalState)
    stream << PhysConst::Instance()->FindParticle(i).GetName() << " ";
  stream << ")";

  LOG(info) << "HelicityKinematics::HelicityKinematics() | Initialize reaction "
            << stream.str();
}

HelicityKinematics::HelicityKinematics(boost::property_tree::ptree pt) {

  auto initialS = pt.get_child("InitialState");
  _initialState = std::vector<int>(initialS.size());
  for (auto i : initialS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = PhysConst::Instance()->FindParticle(name);
    _initialState.at(i.second.get<int>("<xmlattr>.Id")) = partP.GetId();
  }

  auto finalS = pt.get_child("FinalState");
  _finalState = std::vector<int>(finalS.size());
  for (auto i : finalS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = PhysConst::Instance()->FindParticle(name);
    _finalState.at(i.second.get<int>("<xmlattr>.Id")) = partP.GetId();
  }

  auto phspVal = pt.get_optional<double>("PhspVolume");
  if (phspVal) {
    SetPhspVolume(phspVal.get());
  }

  // Creating unique title
  std::stringstream stream;
  stream << "( ";
  for (auto i : _initialState)
    stream << PhysConst::Instance()->FindParticle(i).GetName() << " ";
  stream << ")->( ";
  for (auto i : _finalState)
    stream << PhysConst::Instance()->FindParticle(i).GetName() << " ";
  stream << ")";

  LOG(info) << "HelicityKinematics::HelicityKinematics() | Initialize reaction "
            << stream.str();

  assert(_initialState.size() == 1);
  pid idMother = _initialState.at(0);
  auto motherProp = PhysConst::Instance()->FindParticle(idMother);
  _M = motherProp.GetMass();
  _Msq = _M * _M;
  _spinM = motherProp.GetSpinQuantumNumber("Spin");
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

const std::pair<double, double> &
HelicityKinematics::GetInvMassBounds(const SubSystem sys) const {
  return GetInvMassBounds(
      const_cast<HelicityKinematics *>(this)->GetDataID(sys));
}

const std::pair<double, double> &
HelicityKinematics::GetInvMassBounds(int sysID) const {
  return _invMassBounds.at(sysID);
}

std::pair<double, double>
HelicityKinematics::CalculateInvMassBounds(const SubSystem sys) const {

  /* We use the formulae from (PDG2016 Kinematics Fig.47.3). I hope the
   * generalization to n-body decays is correct.
   */
  std::pair<double, double> lim(0, _M);
  // Sum up masses of all final state particles
  for (auto j : sys.GetFinalStates())
    for (auto i : j)
      lim.first +=
          PhysConst::Instance()->FindParticle(_finalState.at(i)).GetMass();
  lim.first *= lim.first;

  for (auto i : sys.GetRecoilState())
    lim.second -=
        PhysConst::Instance()->FindParticle(_finalState.at(i)).GetMass();
  lim.second *= lim.second;

  return lim;
}

void HelicityKinematics::EventToDataPoint(const Event &event,
                                          dataPoint &point) const {
  assert(_listSubSystem.size() == _invMassBounds.size());

  for (int i = 0; i < _listSubSystem.size(); i++)
    EventToDataPoint(event, point, _listSubSystem.at(i), _invMassBounds.at(i));
  return;
}

void HelicityKinematics::EventToDataPoint(const Event &event, dataPoint &point,
                                          const SubSystem sys) const {
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
    const Event &event, dataPoint &point, const SubSystem sys,
    const std::pair<double, double> limits) const {

  if (sys.GetFinalStates().size() != 2)
    return;

  FourMomentum recoilP4;
  for (auto s : sys.GetRecoilState())
    recoilP4 += event.GetParticle(s).GetFourMomentum();

  FourMomentum finalA, finalB;
  for (auto s : sys.GetFinalStates().at(0))
    finalA += event.GetParticle(s).GetFourMomentum();

  for (auto s : sys.GetFinalStates().at(1))
    finalB += event.GetParticle(s).GetFourMomentum();

  FourMomentum totalP4 = finalA + finalB;
  double mSq = totalP4.GetInvMassSq();

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

  // When using finalB here the WignerD changes sign. In the end this does not
  // matter
  QFT::Vector4<double> qftTotalP4(totalP4);
  QFT::Vector4<double> qftFinalA(finalA);
  QFT::Vector4<double> qftRecoilP4(recoilP4);

  /* Boost one final state four momentum and the four momentum of the recoil
   * system to the center of mass system of the two-body decay
   */
  qftFinalA.Boost(qftTotalP4);
  qftRecoilP4.Boost(qftTotalP4);
  //    qftRecoilP4 *= (-1);

  // Calculate the angles between recoil system and final state.
  qftFinalA.Rotate(qftRecoilP4.Phi(), qftRecoilP4.Theta(),
                   (-1) * qftRecoilP4.Phi());
  double cosTheta = qftFinalA.CosTheta();
  double phi = qftFinalA.Phi();

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
  //sys.GetFinalStates().at(0).at(0)
  //              << sys.GetFinalStates().at(1).at(0) << ") angle ("<<
  //sys.GetFinalStates().at(0).at(0)
  //              << sys.GetRecoilState().at(0) << ") - "
  //              << " (ma=" << ma<<" mb="<<mb<<" mSpec="<<mspec<<"
  //mABSq="<<invMassSqA << " mASpecSq=" << invMassSqB << ") = " << cc
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

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
