// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "qft++/Vector4.h"

#include "Core/DataPoint.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Properties.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

ComPWA::KinematicsProperties
createKinematicsProperties(std::shared_ptr<PartList> partL,
                           const boost::property_tree::ptree &pt) {
  ComPWA::KinematicsProperties KinProps;
  KinProps.ParticleList = partL;
  auto initialS = pt.get_child("InitialState");
  KinProps.InitialState = std::vector<int>(initialS.size());
  for (auto i : initialS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = partL->find(name)->second;
    unsigned int pos = i.second.get<unsigned int>("<xmlattr>.PositionIndex");
    KinProps.InitialState.at(pos) = partP.GetId();
  }
  if (KinProps.InitialState.size() == 1) {
    double sqrtS = FindParticle(partL, KinProps.InitialState.at(0)).GetMass();
    KinProps.InitialStateP4 = ComPWA::FourMomentum(0, 0, 0, sqrtS);
  } else {
    // If the initial state has more than one particle we require that the
    // initial four momentum is specified
    KinProps.InitialStateP4 =
        FourMomentumFactory(pt.get_child("InitialFourMomentum"));
  }

  auto finalS = pt.get_child("FinalState");
  KinProps.FinalState = std::vector<int>(finalS.size());
  KinProps.FinalStateEventPositionMapping =
      std::vector<unsigned int>(finalS.size());
  for (auto i : finalS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = partL->find(name)->second;
    unsigned int id = i.second.get<unsigned int>("<xmlattr>.Id");
    unsigned int pos(id);
    boost::optional<unsigned int> opt_pos =
        i.second.get_optional<unsigned int>("<xmlattr>.PositionIndex");
    if (opt_pos)
      pos = opt_pos.get();
    KinProps.FinalState.at(pos) = partP.GetId();
    KinProps.FinalStateEventPositionMapping[pos] = id;
  }

  return KinProps;
}

HelicityKinematics::HelicityKinematics(std::shared_ptr<PartList> partL,
                                       const std::vector<pid> &initialState,
                                       const std::vector<pid> &finalState,
                                       const ComPWA::FourMomentum &cmsP4)
    : ComPWA::Kinematics(ComPWA::KinematicsProperties(
          initialState, finalState, partL, cmsP4, [&finalState]() {
            std::vector<unsigned int> FinalStateEventPositionMapping;
            for (unsigned int i = 0; i < finalState.size(); ++i) {
              FinalStateEventPositionMapping.push_back(i);
            }
            return FinalStateEventPositionMapping;
          }())) {
  // Create title
  std::stringstream stream;
  const ComPWA::KinematicsProperties &KinProps(getKinematicsProperties());
  stream << "( ";
  for (auto i : KinProps.InitialState)
    stream << FindParticle(partL, i).name() << " ";
  stream << ")->( ";
  for (auto i : KinProps.FinalState)
    stream << FindParticle(partL, i).name() << " ";
  stream << ")";

  LOG(INFO) << "HelicityKinematics::HelicityKinematics() | Initialize reaction "
            << stream.str();
}

HelicityKinematics::HelicityKinematics(std::shared_ptr<PartList> partL,
                                       const boost::property_tree::ptree &pt)
    : Kinematics(createKinematicsProperties(partL, pt)) {
  auto phspVal = pt.get_optional<double>("PhspVolume");
  if (phspVal) {
    setPhspVolume(phspVal.get());
  }
  // Create title
  std::stringstream stream;
  const ComPWA::KinematicsProperties &KinProps(getKinematicsProperties());
  stream << "( ";
  for (auto i : KinProps.InitialState)
    stream << FindParticle(partL, i).name() << " ";
  stream << ")->( ";
  for (auto i : KinProps.FinalState)
    stream << FindParticle(partL, i).name() << " ";
  stream << ")";

  LOG(INFO) << "HelicityKinematics::HelicityKinematics() | Initialize reaction "
            << stream.str();
}

bool HelicityKinematics::isWithinPhsp(const DataPoint &point) const {
  unsigned int subSystemID = 0;
  unsigned int pos = 0;
  while ((pos + 2) < point.size()) {
    auto bounds = invMassBounds(subSystem(subSystemID));
    if (point.value(pos) < bounds.first || point.value(pos) > bounds.second)
      return false;
    if (point.value(pos + 1) < -1 || point.value(pos + 1) > +1)
      return false;
    if (point.value(pos + 2) < 0 || point.value(pos + 2) > 2 * M_PI)
      return false;

    pos += 3;
    subSystemID++;
  }

  return true;
}

void HelicityKinematics::convert(const Event &event, DataPoint &point) const {
  assert(Subsystems.size() == InvMassBounds.size());

  if (!Subsystems.size()) {
    LOG(ERROR) << "HelicityKinematics::convert() | No variables were "
                  "requested before. Therefore this function is doing nothing!";
  }
  for (unsigned int i = 0; i < Subsystems.size(); i++)
    convert(event, point, Subsystems.at(i), InvMassBounds.at(i));
  return;
}

void HelicityKinematics::convert(const Event &event, DataPoint &point,
                                 const SubSystem &sys) const {
  auto massLimits = invMassBounds(sys);
  convert(event, point, sys, massLimits);
}

unsigned int HelicityKinematics::getDataID(const SubSystem &subSys) const {
  auto const result = std::find(Subsystems.begin(), Subsystems.end(), subSys);
  return result - Subsystems.begin();
}

unsigned int HelicityKinematics::addSubSystem(const SubSystem &subSys) {
  // We calculate the variables currently for two-body decays
  if (subSys.getFinalStates().size() != 2) {
    std::stringstream ss;
    ss << "HelicityKinematics::addSubSystem(const SubSystem "
          "&subSys): Number of final state particles = "
       << subSys.getFinalStates().size() << ", which is != 2";
    throw std::runtime_error(ss.str());
  }
  unsigned int pos(Subsystems.size());
  auto const result = std::find(Subsystems.begin(), Subsystems.end(), subSys);
  if (result == Subsystems.end()) {
    Subsystems.push_back(subSys);
    InvMassBounds.push_back(calculateInvMassBounds(subSys));
    std::string subsys_string(subSys.to_string());
    VariableNames.push_back("mSq" + subsys_string);
    VariableNames.push_back("theta" + subsys_string);
    VariableNames.push_back("phi" + subsys_string);
  } else {
    pos = result - Subsystems.begin();
  }
  return pos;
}

unsigned int HelicityKinematics::addSubSystem(
    const std::vector<unsigned int> &FinalA,
    const std::vector<unsigned int> &FinalB,
    const std::vector<unsigned int> &Recoil,
    const std::vector<unsigned int> &ParentRecoil) {
  std::vector<std::vector<unsigned int>> ConvertedFinalStates;
  std::vector<unsigned int> temp;
  const auto &FinalStateEventPositionMapping(
      getKinematicsProperties().FinalStateEventPositionMapping);
  for (auto i : FinalA)
    temp.push_back(FinalStateEventPositionMapping[i]);
  ConvertedFinalStates.push_back(temp);
  temp.clear();
  for (auto i : FinalB)
    temp.push_back(FinalStateEventPositionMapping[i]);
  ConvertedFinalStates.push_back(temp);
  std::vector<unsigned int> ConvertedRecoil;
  for (auto i : Recoil)
    ConvertedRecoil.push_back(FinalStateEventPositionMapping[i]);
  std::vector<unsigned int> ConvertedParentRecoil;
  for (auto i : ParentRecoil)
    ConvertedParentRecoil.push_back(FinalStateEventPositionMapping[i]);

  return addSubSystem(
      SubSystem(ConvertedFinalStates, ConvertedRecoil, ConvertedParentRecoil));
}

double HelicityKinematics::helicityAngle(double M, double m, double m2,
                                         double mSpec, double invMassSqA,
                                         double invMassSqB) const {
  // Calculate energy and momentum of m1/m2 in the invMassSqA rest frame
  double eCms = (invMassSqA + m * m - m2 * m2) / (2.0 * std::sqrt(invMassSqA));
  double pCms = eCms * eCms - m * m;
  // Calculate energy and momentum of mSpec in the invMassSqA rest frame
  double eSpecCms =
      (M * M - mSpec * mSpec - invMassSqA) / (2.0 * std::sqrt(invMassSqA));
  double pSpecCms = eSpecCms * eSpecCms - mSpec * mSpec;
  double cosAngle =
      -(invMassSqB - m * m - mSpec * mSpec - 2.0 * eCms * eSpecCms) /
      (2 * std::sqrt(pCms * pSpecCms));

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

void HelicityKinematics::convert(const Event &event, DataPoint &point,
                                 const SubSystem &sys,
                                 const std::pair<double, double> limits) const {
  FourMomentum FinalA, FinalB;
  for (auto s : sys.getFinalStates().at(0)) {
    unsigned int index = convertFinalStateIDToPositionIndex(s);
    FinalA += event.particle(index).fourMomentum();
  }

  for (auto s : sys.getFinalStates().at(1)) {
    unsigned int index = convertFinalStateIDToPositionIndex(s);
    FinalB += event.particle(index).fourMomentum();
  }

  // Four momentum of the decaying resonance
  FourMomentum State = FinalA + FinalB;
  double mSq = State.invMassSq();

  if (mSq <= limits.first) {
    // We allow for a deviation from the limits of 10 times the numerical
    // precision
    if (ComPWA::equal(mSq, limits.first, 10))
      mSq = limits.first;
    else
      throw BeyondPhsp("HelicityKinematics::convert() |"
                       " Point beyond phase space boundaries!");
  }
  if (mSq >= limits.second) {
    // We allow for a deviation from the limits of 10 times the numerical
    // precision
    if (ComPWA::equal(mSq, limits.second, 10))
      mSq = limits.second;
    else
      throw BeyondPhsp("HelicityKinematics::convert() |"
                       " Point beyond phase space boundaries!");
  }

  QFT::Vector4<double> DecayingState(State);
  QFT::Vector4<double> Daughter(FinalA);

  // the first step is boosting everything into the rest system of the
  // decaying state
  Daughter.Boost(DecayingState);

  // calculate the recoil and parent recoil
  auto const &RecoilState = sys.getRecoilState();
  if (RecoilState.size() > 0) {
    FourMomentum TempRecoil;
    for (auto s : RecoilState) {
      unsigned int index = convertFinalStateIDToPositionIndex(s);
      TempRecoil += event.particle(index).fourMomentum();
    }
    QFT::Vector4<double> Recoil(TempRecoil);
    Recoil.Boost(DecayingState);

    // rotate vectors so that recoil moves in the negative z-axis direction
    Daughter.RotateZ(-Recoil.Phi());
    Daughter.RotateY(M_PI - Recoil.Theta());

    auto const &ParentRecoilState = sys.getParentRecoilState();
    QFT::Vector4<double> ParentRecoil;
    if (ParentRecoilState.size() > 0) {
      FourMomentum TempParentRecoil;
      for (auto s : ParentRecoilState) {
        unsigned int index = convertFinalStateIDToPositionIndex(s);
        TempParentRecoil += event.particle(index).fourMomentum();
      }
      ParentRecoil = TempParentRecoil;
    } else {
      // in case there is no parent recoil, it is artificially along z
      ParentRecoil.SetP4(0, 0, 0, 1.0);
    }

    ParentRecoil.Boost(DecayingState);

    ParentRecoil.RotateZ(-Recoil.Phi());
    ParentRecoil.RotateY(M_PI - Recoil.Theta());

    // rotate around the z-axis so that the parent recoil lies in the x-z
    // plane
    Daughter.RotateZ(M_PI - ParentRecoil.Phi());
  }

  double cosTheta = Daughter.CosTheta();
  double phi = Daughter.Phi();

  //   Check if values are within allowed range.
  if (cosTheta > 1 || cosTheta < -1 || phi > M_PI || phi < (-1) * M_PI ||
      std::isnan(cosTheta) || std::isnan(phi)) {
    throw BeyondPhsp("HelicityKinematics::convert() |"
                     " Point beypond phase space boundaries!");
  }

  point.setWeight(event.weight());
  point.setEfficiency(event.efficiency());
  point.values().push_back(mSq);
  point.values().push_back(std::acos(cosTheta));
  point.values().push_back(phi);
}

const std::pair<double, double> &
HelicityKinematics::invMassBounds(const SubSystem &sys) const {
  return invMassBounds(getDataID(sys));
}

const std::pair<double, double> &
HelicityKinematics::invMassBounds(int sysID) const {
  return InvMassBounds.at(sysID);
}

std::pair<double, double>
HelicityKinematics::calculateInvMassBounds(const SubSystem &sys) const {
  /// We use the formulae from (PDG2016 Kinematics Fig.47.3). I hope the
  /// generalization to n-body decays is correct.
  const auto &KinProps(getKinematicsProperties());
  std::pair<double, double> lim(0, KinProps.InitialStateP4.invMass());
  // Sum up masses of all final state particles
  for (auto j : sys.getFinalStates())
    for (auto i : j) {
      unsigned int index = convertFinalStateIDToPositionIndex(i);
      lim.first += FindParticle(KinProps.ParticleList, KinProps.FinalState.at(index)).GetMass();
    }
  lim.first *= lim.first;

  for (auto i : sys.getRecoilState()) {
    unsigned int index = convertFinalStateIDToPositionIndex(i);
    lim.second -= FindParticle(KinProps.ParticleList, KinProps.FinalState.at(index)).GetMass();
  }
  lim.second *= lim.second;

  return lim;
}

} // ns::HelicityFormalism
} // ns::Physics
} // ns::ComPWA
