// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <algorithm>
#include <cmath>
#include <numeric>

#include "Core/DataPoint.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Properties.hpp"

#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/qft++/Vector4.h"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

HelicityKinematics::HelicityKinematics(std::shared_ptr<PartList> partL,
                                       std::vector<pid> initialState,
                                       std::vector<pid> finalState,
                                       ComPWA::FourMomentum cmsP4)
    : Kinematics(initialState, finalState, cmsP4), ParticleList(partL) {
  FinalStateEventPositionMapping.clear();
  for (unsigned int i = 0; i < FinalState.size(); ++i) {
    FinalStateEventPositionMapping.push_back(i);
  }

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
  for (auto i : FinalState)
    stream << FindParticle(partL, i).name() << " ";
  stream << ")";

  LOG(info) << "HelicityKinematics::HelicityKinematics() | Initialize reaction "
            << stream.str();
  return;
}

HelicityKinematics::HelicityKinematics(std::shared_ptr<PartList> partL,
                                       boost::property_tree::ptree pt)
    : ParticleList(partL) {
  auto initialS = pt.get_child("InitialState");
  InitialState = std::vector<int>(initialS.size());
  for (auto i : initialS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = partL->find(name)->second;
    unsigned int pos = i.second.get<unsigned int>("<xmlattr>.PositionIndex");
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
  FinalStateEventPositionMapping = std::vector<unsigned int>(finalS.size());
  for (auto i : finalS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = partL->find(name)->second;
    unsigned int id = i.second.get<unsigned int>("<xmlattr>.Id");
    unsigned int pos(id);
    boost::optional<unsigned int> opt_pos =
        i.second.get_optional<unsigned int>("<xmlattr>.PositionIndex");
    if (opt_pos)
      pos = opt_pos.get();
    FinalState.at(pos) = partP.GetId();
    FinalStateEventPositionMapping[pos] = id;
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

  LOG(info) << "HelicityKinematics::HelicityKinematics() | Initialize reaction "
            << stream.str();
  return;
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
    LOG(error) << "HelicityKinematics::convert() | No variabels were "
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

unsigned int HelicityKinematics::convertFinalStateIDToPositionIndex(
    unsigned int fs_id) const {
  auto result = std::find(FinalStateEventPositionMapping.begin(),
                          FinalStateEventPositionMapping.end(), fs_id);
  return std::distance(FinalStateEventPositionMapping.begin(), result);
}

std::vector<unsigned int>
HelicityKinematics::convertFinalStateIDToPositionIndex(
    const std::vector<unsigned int> &fs_ids) const {
  std::vector<unsigned int> pos_indices;
  pos_indices.reserve(fs_ids.size());
  for (auto fs_id : fs_ids) {
    pos_indices.push_back(convertFinalStateIDToPositionIndex(fs_id));
  }
  return pos_indices;
}

unsigned int HelicityKinematics::getDataID(const SubSystem &subSys) const {
  auto const result = std::find(Subsystems.begin(), Subsystems.end(), subSys);
  return result - Subsystems.begin();
}

unsigned int HelicityKinematics::addSubSystem(const SubSystem &subSys) {
  // We calculate the variables currently for two-body decays
  if (subSys.GetFinalStates().size() != 2)
    return 0;
  unsigned int pos(Subsystems.size());
  auto const result = std::find(Subsystems.begin(), Subsystems.end(), subSys);
  if (result == Subsystems.end()) {
    Subsystems.push_back(subSys);
    InvMassBounds.push_back(calculateInvMassBounds(subSys));
    std::string subsys_string(subSys.to_string());
    VariableNames.push_back("mSq_" + subsys_string);
    VariableNames.push_back("cosTheta_" + subsys_string);
    VariableNames.push_back("phi_" + subsys_string);
    VariableTitles.push_back("m^{2}(" + subsys_string + ")");
    VariableTitles.push_back("#cos#Theta(" + subsys_string + ")");
    VariableTitles.push_back("#phi(" + subsys_string + ")");
  } else {
    pos = result - Subsystems.begin();
  }
  return pos;
}

unsigned int
HelicityKinematics::addSubSystem(const std::vector<unsigned int> &recoilS,
                                 const std::vector<unsigned int> &finalA,
                                 const std::vector<unsigned int> &finalB) {
  std::vector<unsigned int> conv_recoilS;
  std::vector<std::vector<unsigned int>> conv_finalS;
  for (auto i : recoilS)
    conv_recoilS.push_back(FinalStateEventPositionMapping[i]);
  std::vector<unsigned int> temp;
  for (auto i : finalA)
    temp.push_back(FinalStateEventPositionMapping[i]);
  conv_finalS.push_back(temp);
  temp.clear();
  for (auto i : finalB)
    temp.push_back(FinalStateEventPositionMapping[i]);
  conv_finalS.push_back(temp);

  return addSubSystem(SubSystem(conv_recoilS, conv_finalS));
}

double HelicityKinematics::helicityAngle(double M, double m, double m2,
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

void HelicityKinematics::convert(const Event &event, DataPoint &point,
                                 const SubSystem &sys,
                                 const std::pair<double, double> limits) const {
  assert(sys.GetFinalStates().size() == 2 &&
         "HelicityKinematics::convert() | More then two particles.");

  FourMomentum cms;
  for (auto s : sys.GetRecoilState()) {
    unsigned int index = convertFinalStateIDToPositionIndex(s);
    cms += event.particle(index).fourMomentum();
  }

  FourMomentum finalA, finalB;
  for (auto s : sys.GetFinalStates().at(0)) {
    unsigned int index = convertFinalStateIDToPositionIndex(s);
    finalA += event.particle(index).fourMomentum();
  }

  for (auto s : sys.GetFinalStates().at(1)) {
    unsigned int index = convertFinalStateIDToPositionIndex(s);
    finalB += event.particle(index).fourMomentum();
  }

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
  }

  point.values().push_back(mSq);
  point.values().push_back(cosTheta);
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
  std::pair<double, double> lim(0, InitialStateP4.invMass());
  // Sum up masses of all final state particles
  for (auto j : sys.GetFinalStates())
    for (auto i : j) {
      unsigned int index = convertFinalStateIDToPositionIndex(i);
      lim.first += FindParticle(ParticleList, FinalState.at(index)).GetMass();
    }
  lim.first *= lim.first;

  for (auto i : sys.GetRecoilState()) {
    unsigned int index = convertFinalStateIDToPositionIndex(i);
    lim.second -= FindParticle(ParticleList, FinalState.at(index)).GetMass();
  }
  lim.second *= lim.second;

  return lim;
}

} // ns::HelicityFormalism
} // ns::Physics
} // ns::ComPWA
