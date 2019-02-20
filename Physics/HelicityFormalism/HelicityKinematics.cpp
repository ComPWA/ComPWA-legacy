// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <tuple>

#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Properties.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

#include "qft++/Vector4.h"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

HelicityKinematics::HelicityKinematics(std::shared_ptr<PartList> partL,
                                       const std::vector<pid> &initialState,
                                       const std::vector<pid> &finalState,
                                       const ComPWA::FourMomentum &cmsP4)
    : HelicityKinematics(ParticleStateTransitionKinematicsInfo(
          initialState, finalState, partL, cmsP4, [&finalState]() {
            std::vector<unsigned int> FinalStateEventPositionMapping;
            for (unsigned int i = 0; i < finalState.size(); ++i) {
              FinalStateEventPositionMapping.push_back(i);
            }
            return FinalStateEventPositionMapping;
          }())) {}

HelicityKinematics::HelicityKinematics(
    const ParticleStateTransitionKinematicsInfo &ki)
    : HelicityKinematics(ki, 1.0) {
  ///  Calculation of n-dimensional phase space volume.
  ///  ToDo: We need to implement an analytical calculation here
}

HelicityKinematics::HelicityKinematics(
    const ParticleStateTransitionKinematicsInfo &ki, double phspvol)
    : KinematicsInfo(ki), PhspVolume(phspvol) {
  LOG(INFO) << "HelicityKinematics::"
               "HelicityKinematics() | Initialized kinematics "
               "for reaction "
            << KinematicsInfo;
}

double HelicityKinematics::phspVolume() const { return PhspVolume; }

bool HelicityKinematics::isWithinPhaseSpace(const DataPoint &point) const {
  unsigned int subSystemID = 0;
  unsigned int pos = 0;
  while ((pos + 2) < point.KinematicVariableList.size()) {
    auto bounds = invMassBounds(subSystem(subSystemID));
    if (point.KinematicVariableList[pos] < bounds.first ||
        point.KinematicVariableList[pos] > bounds.second)
      return false;
    if (point.KinematicVariableList[pos + 1] < 0 ||
        point.KinematicVariableList[pos + 1] > M_PI)
      return false;
    if (point.KinematicVariableList[pos + 2] < -M_PI ||
        point.KinematicVariableList[pos + 2] > M_PI)
      return false;

    pos += 3;
    subSystemID++;
  }

  return true;
}

DataPoint HelicityKinematics::convert(const Event &event) const {
  assert(Subsystems.size() == InvMassBounds.size());

  if (!Subsystems.size()) {
    LOG(ERROR) << "HelicityKinematics::convert() | No variables were "
                  "requested before. Therefore this function is doing nothing!";
  }
  DataPoint point;
  for (unsigned int i = 0; i < Subsystems.size(); i++)
    convert(event, point, Subsystems.at(i), InvMassBounds.at(i));
  return point;
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

void HelicityKinematics::createAllSubsystems() {
  std::vector<IndexListTuple> CurrentSubsystems;
  LOG(INFO) << "creating all Subsystems!";

  std::vector<IndexListTuple> AllSubsystems;
  // add current subsystems
  for (auto const &SubSys : subSystems()) {
    AllSubsystems.push_back(
        std::make_tuple(KinematicsInfo.convertFinalStateIDToPositionIndex(
                            SubSys.getFinalStates()[0]),
                        KinematicsInfo.convertFinalStateIDToPositionIndex(
                            SubSys.getFinalStates()[1]),
                        KinematicsInfo.convertFinalStateIDToPositionIndex(
                            SubSys.getRecoilState()),
                        KinematicsInfo.convertFinalStateIDToPositionIndex(
                            SubSys.getParentRecoilState())));
  }

  IndexList FinalStateIDs(KinematicsInfo.getFinalStateParticleCount());
  unsigned int i = 0;
  std::generate(FinalStateIDs.begin(), FinalStateIDs.end(),
                [&i]() -> unsigned int { return i++; });
  for (auto const &x : redistributeIndexLists(FinalStateIDs, IndexList())) {
    CurrentSubsystems.push_back(
        std::make_tuple(x.first, x.second, IndexList(), IndexList()));
  }

  while (CurrentSubsystems.size() > 0) {
    auto TempSubsystems(CurrentSubsystems);
    CurrentSubsystems.clear();
    for (auto const &x : TempSubsystems) {
      IndexList CurrentA(std::get<0>(x));
      IndexList CurrentB(std::get<1>(x));
      IndexList CurrentRecoil(std::get<2>(x));
      IndexList CurrentParentRecoil(std::get<3>(x));
      auto SortedEntry = sortSubsystem(x);
      auto result = std::find_if(
          AllSubsystems.begin(), AllSubsystems.end(),
          [&SortedEntry, this](const IndexListTuple &element) -> bool {
            return sortSubsystem(element) == SortedEntry;
          });
      if (result == AllSubsystems.end()) {
        AllSubsystems.push_back(SortedEntry);
      }

      // create more combinations on this level
      for (auto const &ABnew : redistributeIndexLists(CurrentA, CurrentB)) {
        CurrentSubsystems.push_back(std::make_tuple(
            ABnew.first, ABnew.second, CurrentRecoil, CurrentParentRecoil));
      }
      // Try to go a level deeper for A
      if (CurrentA.size() > 1) {
        for (auto const &ABnew :
             redistributeIndexLists(CurrentA, IndexList())) {
          // use current B as recoil
          // and shift recoil to parent recoil
          CurrentSubsystems.push_back(std::make_tuple(ABnew.first, ABnew.second,
                                                      CurrentB, CurrentRecoil));
        }
      }
      // same for B
      if (CurrentB.size() > 1) {
        for (auto const &ABnew :
             redistributeIndexLists(CurrentB, IndexList())) {
          CurrentSubsystems.push_back(std::make_tuple(ABnew.first, ABnew.second,
                                                      CurrentA, CurrentRecoil));
        }
      }
    }
  }
  for (auto const &x : AllSubsystems) {
    addSubSystem(std::get<0>(x), std::get<1>(x), std::get<2>(x),
                 std::get<3>(x));
  }
}

HelicityKinematics::IndexListTuple
HelicityKinematics::sortSubsystem(const IndexListTuple &SubSys) const {
  IndexList FinalStateA(std::get<0>(SubSys));
  IndexList FinalStateB(std::get<1>(SubSys));
  IndexList RecoilState(std::get<2>(SubSys));
  IndexList ParentRecoilState(std::get<3>(SubSys));

  // create sorted entry
  std::sort(FinalStateA.begin(), FinalStateA.end());
  std::sort(FinalStateB.begin(), FinalStateB.end());
  std::sort(RecoilState.begin(), RecoilState.end());
  std::sort(ParentRecoilState.begin(), ParentRecoilState.end());

  IndexListTuple SortedTuple;
  if (FinalStateA > FinalStateB)
    SortedTuple = std::make_tuple(FinalStateB, FinalStateA, RecoilState,
                                  ParentRecoilState);
  else
    SortedTuple = std::make_tuple(FinalStateA, FinalStateB, RecoilState,
                                  ParentRecoilState);
  return SortedTuple;
}

std::vector<std::pair<IndexList, IndexList>>
HelicityKinematics::redistributeIndexLists(const IndexList &A,
                                           const IndexList &B) const {
  std::vector<std::pair<IndexList, IndexList>> NewIndexLists;
  if (A.size() < B.size())
    throw std::runtime_error("HelicityKinematics::redistributeIndexLists(): A "
                             "cannot have less content than B!");
  if (A.size() == 2 && B.size() == 0) {
    NewIndexLists.push_back(std::make_pair(IndexList{A[0]}, IndexList{A[1]}));
  } else if (A.size() - B.size() > 1) {
    for (unsigned int i = 0; i < A.size(); ++i) {
      IndexList TempB(B);
      TempB.push_back(A[i]);
      IndexList TempA(A);
      TempA.erase(TempA.begin() + i);
      NewIndexLists.push_back(std::make_pair(TempA, TempB));
    }
  }
  return NewIndexLists;
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
    std::stringstream ss;
    ss << subSys;
    VariableNames.push_back("mSq" + ss.str());
    VariableNames.push_back("theta" + ss.str());
    VariableNames.push_back("phi" + ss.str());
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
  ConvertedFinalStates.push_back(
      KinematicsInfo.convertPositionIndexToFinalStateID(FinalA));

  ConvertedFinalStates.push_back(
      KinematicsInfo.convertPositionIndexToFinalStateID(FinalB));

  std::vector<unsigned int> ConvertedRecoil =
      KinematicsInfo.convertPositionIndexToFinalStateID(Recoil);
  std::vector<unsigned int> ConvertedParentRecoil =
      KinematicsInfo.convertPositionIndexToFinalStateID(ParentRecoil);

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
      (2.0 * std::sqrt(pCms * pSpecCms));
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
  /*if (cosAngle > 1.0 || cosAngle < -1.0) { // faster
    if (ComPWA::equal(cosAngle, 1.0, 10))
      cosAngle = 1.0;
    else if (ComPWA::equal(cosAngle, -1.0, 10))
      cosAngle = -1.0;
    else {
      throw BeyondPhsp("HelicityKinematics::helicityAngle() | "
                       "scattering angle out of range! Datapoint beyond"
                       "phsp?");
    }
  }*/
  return cosAngle;
}

void HelicityKinematics::convert(const Event &event, DataPoint &point,
                                 const SubSystem &sys,
                                 const std::pair<double, double> limits) const {
  FourMomentum FinalA, FinalB;
  for (auto s : sys.getFinalStates().at(0)) {
    unsigned int index = KinematicsInfo.convertFinalStateIDToPositionIndex(s);
    FinalA += event.ParticleList[index].fourMomentum();
  }

  for (auto s : sys.getFinalStates().at(1)) {
    unsigned int index = KinematicsInfo.convertFinalStateIDToPositionIndex(s);
    FinalB += event.ParticleList[index].fourMomentum();
  }

  // Four momentum of the decaying resonance
  FourMomentum State = FinalA + FinalB;
  double mSq = State.invMassSq();

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
      unsigned int index = KinematicsInfo.convertFinalStateIDToPositionIndex(s);
      TempRecoil += event.ParticleList[index].fourMomentum();
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
        unsigned int index =
            KinematicsInfo.convertFinalStateIDToPositionIndex(s);
        TempParentRecoil += event.ParticleList[index].fourMomentum();
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

  point.Weight = event.Weight;
  point.KinematicVariableList.push_back(mSq);
  point.KinematicVariableList.push_back(std::acos(cosTheta));
  point.KinematicVariableList.push_back(phi);
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
  std::pair<double, double> lim(0.0,
                                KinematicsInfo.getInitialStateInvariantMass());
  // Sum up masses of all final state particles
  for (auto j : sys.getFinalStates())
    lim.first += KinematicsInfo.calculateFinalStateIDMassSum(j);
  lim.first *= lim.first;
  // we add a space of numeric double precision to the boundary
  lim.first -= std::numeric_limits<double>::epsilon() * lim.first;
  lim.second -=
      KinematicsInfo.calculateFinalStateIDMassSum(sys.getRecoilState());
  lim.second *= lim.second;
  lim.second += std::numeric_limits<double>::epsilon() * lim.second;

  return lim;
}

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA
