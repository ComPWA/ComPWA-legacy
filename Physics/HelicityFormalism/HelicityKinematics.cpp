// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <tuple>

#include "Core/Event.hpp"
#include "Core/Exceptions.hpp"
#include "Core/Logging.hpp"
#include "Core/FourMomentum.hpp"
#include "Core/Properties.hpp"
#include "Data/DataSet.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

#include "qft++/Vector4.h"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

HelicityKinematics::HelicityKinematics(ComPWA::ParticleList partL,
                                       std::vector<pid> initialState,
                                       std::vector<pid> finalState,
                                       ComPWA::FourMomentum cmsP4)
    : HelicityKinematics(ParticleStateTransitionKinematicsInfo(
          initialState, finalState, partL, cmsP4, [&finalState]() {
            std::vector<unsigned int> FinalStateEventPositionMapping;
            for (unsigned int i = 0; i < finalState.size(); ++i) {
              FinalStateEventPositionMapping.push_back(i);
            }
            return FinalStateEventPositionMapping;
          }())) {}

HelicityKinematics::HelicityKinematics(ParticleStateTransitionKinematicsInfo ki)
    : HelicityKinematics(ki, 1.0) {
  ///  Calculation of n-dimensional phase space volume.
  ///  ToDo: We need to implement an analytical calculation here
}

HelicityKinematics::HelicityKinematics(ParticleStateTransitionKinematicsInfo ki,
                                       double PhspVol)
    : KinematicsInfo(ki), PhspVolume(PhspVol) {
  LOG(INFO) << "HelicityKinematics::"
               "HelicityKinematics() | Initialized kinematics "
               "for reaction "
            << KinematicsInfo;
}

double HelicityKinematics::phspVolume() const { return PhspVolume; }

EventCollection HelicityKinematics::reduceToPhaseSpace(
    const EventCollection &DataSample) const {
  EventCollection PhspSample{DataSample.Pids};
  LOG(INFO) << "HelicityKinematics::reduceToPhaseSpace(): "
               "Remove all events outside PHSP boundary from data sample.";

  std::copy_if(DataSample.Events.begin(), DataSample.Events.end(),
               std::back_inserter(PhspSample.Events), [this](auto evt) {
                 for (auto const &x : InvariantMassesSquared) {
                   auto bounds = InvMassBounds.at(x.second);
                   double mass =
                       this->calculateInvariantMassSquared(evt, x.first);
                   if (mass < bounds.first || mass > bounds.second)
                     return false;
                 }
                 for (auto const &x : Subsystems) {
                   auto angles = this->calculateHelicityAngles(evt, x.first);
                   if (angles.first < 0 || angles.first > M_PI)
                     return false;
                   if (angles.second < -M_PI || angles.second > M_PI)
                     return false;
                 }
                 return true;
               });

  LOG(INFO) << "reduceToPhaseSpace(): Removed "
            << DataSample.Events.size() - PhspSample.Events.size() << " from "
            << DataSample.Events.size() << " Events ("
            << (1.0 - PhspSample.Events.size() / DataSample.Events.size()) * 100
            << "%).";

  return PhspSample;
}

std::pair<double, double>
HelicityKinematics::calculateHelicityAngles(const Event &Event,
                                            const SubSystem &SubSys) const {
  FourMomentum FinalA, FinalB;
  for (auto s : SubSys.getFinalStates().at(0)) {
    unsigned int index = KinematicsInfo.convertFinalStateIDToPositionIndex(s);
    FinalA += Event.FourMomenta[index];
  }

  for (auto s : SubSys.getFinalStates().at(1)) {
    unsigned int index = KinematicsInfo.convertFinalStateIDToPositionIndex(s);
    FinalB += Event.FourMomenta[index];
  }

  // Four momentum of the decaying resonance
  FourMomentum State = FinalA + FinalB;

  QFT::Vector4<double> DecayingState(State);
  QFT::Vector4<double> Daughter(FinalA);

  // calculate the recoil and parent recoil
  auto const &RecoilState = SubSys.getRecoilState();

  // the first step is boosting everything into the rest system of the
  // decaying state
  Daughter.Boost(DecayingState);

  if (RecoilState.size() > 0) {
    FourMomentum TempRecoil;
    for (auto s : RecoilState) {
      unsigned int index = KinematicsInfo.convertFinalStateIDToPositionIndex(s);
      TempRecoil += Event.FourMomenta[index];
    }
    QFT::Vector4<double> Recoil(TempRecoil);

    Recoil.Boost(DecayingState);

    // rotate recoil so that recoil points in the -z direction
    Daughter.RotateZ(-Recoil.Phi());
    Daughter.RotateY(M_PI - Recoil.Theta());

    auto const &ParentRecoilState = SubSys.getParentRecoilState();
    // in case there is no parent recoil, it is artificially along z
    QFT::Vector4<double> ParentRecoil(0.0, 0.0, 0.0, 1.0);
    if (ParentRecoilState.size() > 0) {
      FourMomentum TempParentRecoil;
      for (auto s : ParentRecoilState) {
        unsigned int index =
            KinematicsInfo.convertFinalStateIDToPositionIndex(s);
        TempParentRecoil += Event.FourMomenta[index];
      }
      ParentRecoil = TempParentRecoil;
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

  return std::make_pair(std::acos(cosTheta), phi);
}

double HelicityKinematics::calculateInvariantMassSquared(
    const Event &Event, const IndexList &FinalStateIDs) const {
  FourMomentum State;
  for (auto s : FinalStateIDs) {
    unsigned int index = KinematicsInfo.convertFinalStateIDToPositionIndex(s);
    State += Event.FourMomenta[index];
  }

  return State.invariantMassSquared();
}

ComPWA::Data::DataSet
HelicityKinematics::convert(const EventCollection &DataSample) const {
  ComPWA::Data::DataSet Dataset;
  if (!Subsystems.size()) {
    LOG(ERROR) << "HelicityKinematics::convert() | No variables were requested "
                  "before. Therefore this function is doing nothing!";
    return Dataset;
  }
  if (KinematicsInfo.getFinalStatePIDs() != DataSample.Pids) {
    std::stringstream Message;
    Message << "Pids in EventCollection and in Kinematics do not match";
    Message << std::endl << "  ";
    Message << DataSample.Pids.size() << " PIDs in EventCollection:";
    for (auto Pid : DataSample.Pids)
      Message << " " << Pid;
    Message << std::endl << "  ";
    Message << KinematicsInfo.getFinalStatePIDs().size()
            << " PIDs in Kinematics:";
    for (auto Pid : KinematicsInfo.getFinalStatePIDs())
      Message << " " << Pid;
    throw ComPWA::BadParameter(Message.str());
  }
  if (!DataSample.checkPidMatchesEvents()) {
    throw ComPWA::BadParameter("HelicityKinematics::convert() | number of PIDs "
                               "not equal to number of four-momenta");
  }

  for (auto const &Masses : InvariantMassesSquared) {
    Dataset.Data.insert(std::make_pair(Masses.second, std::vector<double>()));
    for (auto const &Event : DataSample.Events) {
      auto Mass = calculateInvariantMassSquared(Event, Masses.first);
      Dataset.Data.at(Masses.second).push_back(Mass);
    }
  }
  for (auto const &SubSystem : Subsystems) {
    Dataset.Data.insert(
        std::make_pair(SubSystem.second.first, std::vector<double>()));
    Dataset.Data.insert(
        std::make_pair(SubSystem.second.second, std::vector<double>()));
    for (auto const &event : DataSample.Events) {
      auto Angles = calculateHelicityAngles(event, SubSystem.first);
      Dataset.Data.at(SubSystem.second.first).push_back(Angles.first);
      Dataset.Data.at(SubSystem.second.second).push_back(Angles.second);
    }
  }
  for (auto const &Event : DataSample.Events) {
    Dataset.Weights.push_back(Event.Weight);
  }
  return Dataset;
}

std::string
HelicityKinematics::registerInvariantMassSquared(IndexList MomentaIDs) {
  std::sort(MomentaIDs.begin(), MomentaIDs.end());
  std::stringstream Name;
  Name << "mSq_(";
  std::string comma("");
  for (auto x : MomentaIDs) {
    Name << comma << x;
    comma = ",";
  }
  Name << ")";

  InvMassBounds.insert(
      std::make_pair(Name.str(), calculateInvMassBounds(MomentaIDs)));
  InvariantMassesSquared.insert(std::make_pair(MomentaIDs, Name.str()));
  return Name.str();
}

std::pair<std::string, std::string>
HelicityKinematics::registerHelicityAngles(SubSystem SubSys) {
  std::stringstream ss;
  ss << SubSys;

  std::string ThetaName("theta" + ss.str());
  std::string PhiName("phi" + ss.str());

  Subsystems.insert(std::make_pair(SubSys, std::make_pair(ThetaName, PhiName)));
  return std::make_pair(ThetaName, PhiName);
}

std::vector<std::pair<ComPWA::IndexList, ComPWA::IndexList>>
redistributeIndexLists(const ComPWA::IndexList &A, const ComPWA::IndexList &B) {
  using ComPWA::IndexList;
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

using IndexListTuple = std::tuple<IndexList, IndexList, IndexList, IndexList>;

IndexListTuple sortSubsystem(const IndexListTuple &SubSys) {
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

void HelicityKinematics::createAllSubsystems() {
  std::vector<IndexListTuple> CurrentSubsystems;
  LOG(INFO) << "creating all Subsystems!";

  std::vector<IndexListTuple> AllSubsystems;
  // add current subsystems
  for (auto const &SubSys : Subsystems) {
    AllSubsystems.push_back(
        std::make_tuple(KinematicsInfo.convertFinalStateIDToPositionIndex(
                            SubSys.first.getFinalStates()[0]),
                        KinematicsInfo.convertFinalStateIDToPositionIndex(
                            SubSys.first.getFinalStates()[1]),
                        KinematicsInfo.convertFinalStateIDToPositionIndex(
                            SubSys.first.getRecoilState()),
                        KinematicsInfo.convertFinalStateIDToPositionIndex(
                            SubSys.first.getParentRecoilState())));
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
      for (auto const &ABNew : redistributeIndexLists(CurrentA, CurrentB)) {
        CurrentSubsystems.push_back(std::make_tuple(
            ABNew.first, ABNew.second, CurrentRecoil, CurrentParentRecoil));
      }
      // Try to go a level deeper for A
      if (CurrentA.size() > 1) {
        for (auto const &ABNew :
             redistributeIndexLists(CurrentA, IndexList())) {
          // use current B as recoil
          // and shift recoil to parent recoil
          CurrentSubsystems.push_back(std::make_tuple(ABNew.first, ABNew.second,
                                                      CurrentB, CurrentRecoil));
        }
      }
      // same for B
      if (CurrentB.size() > 1) {
        for (auto const &ABNew :
             redistributeIndexLists(CurrentB, IndexList())) {
          CurrentSubsystems.push_back(std::make_tuple(ABNew.first, ABNew.second,
                                                      CurrentA, CurrentRecoil));
        }
      }
    }
  }
  for (auto const &x : AllSubsystems) {
    registerSubSystem(std::get<0>(x), std::get<1>(x), std::get<2>(x),
                      std::get<3>(x));
  }
}

std::tuple<std::string, std::string, std::string>
HelicityKinematics::registerSubSystem(const SubSystem &NewSys) {
  // We calculate the variables currently for two-body decays
  if (NewSys.getFinalStates().size() != 2) {
    std::stringstream ss;
    ss << "HelicityKinematics::registerSubSystem(const SubSystem "
          "&subSys): Number of final state particles = "
       << NewSys.getFinalStates().size() << ", which is != 2";
    throw std::runtime_error(ss.str());
  }

  IndexList FS1 = NewSys.getFinalStates().at(0);
  IndexList FS2 = NewSys.getFinalStates().at(1);
  FS1.insert(FS1.end(), FS2.begin(), FS2.end());

  auto MassName = registerInvariantMassSquared(FS1);
  auto AngleNames = registerHelicityAngles(NewSys);

  return std::make_tuple(MassName, AngleNames.first, AngleNames.second);
}

std::tuple<std::string, std::string, std::string>
HelicityKinematics::registerSubSystem(
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

  return registerSubSystem(
      SubSystem(ConvertedFinalStates, ConvertedRecoil, ConvertedParentRecoil));
}

const std::pair<double, double> &HelicityKinematics::getInvariantMassBounds(
    const std::string &InvariantMassName) const {
  return InvMassBounds.at(InvariantMassName);
}

std::pair<double, double> HelicityKinematics::calculateInvMassBounds(
    const IndexList &FinalStateIDs) const {
  /// We use the formulae from (PDG2016 Kinematics Fig.47.3). I hope the
  /// generalization to n-body decays is correct.
  std::pair<double, double> lim(0.0, 0.0);
  // Sum up masses of all final state particles
  lim.first = KinematicsInfo.calculateFinalStateIDMassSum(FinalStateIDs);

  double S = KinematicsInfo.getInitialStateInvariantMassSquared();
  double RemainderMass(0.0);
  for (auto x : KinematicsInfo.getFinalStateMasses())
    RemainderMass += x;
  RemainderMass -= lim.first;

  lim.first *= lim.first;
  // to improve precision: (M - m)^2 -> S - 2 sqrt(S) m + m^2 with S=M^2
  lim.second =
      S - 2 * std::sqrt(S) * RemainderMass + RemainderMass * RemainderMass;

  // extend the invariant mass interval by the numeric double precision
  // otherwise quite a few events at the phase space boundary can be lost
  // due to numerical imprecision of event generators
  // (especially drastic effect, if gammas are in the final state)
  lim.first -= std::numeric_limits<double>::epsilon() * lim.first;
  lim.second += std::numeric_limits<double>::epsilon() * lim.second;
  return lim;
}

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA
