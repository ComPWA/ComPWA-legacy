// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Physics/EvtGen/DalitzKinematics.hpp"

#include "Core/Event.hpp"
#include "Core/Logging.hpp"
#include "Core/FourMomentum.hpp"
#include "Data/DataSet.hpp"

#include "ThirdParty/qft++/include/qft++/Vector4.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>

namespace ComPWA {
namespace Physics {
namespace EvtGen {

DalitzKinematics::DalitzKinematics(ComPWA::ParticleList partL,
                                   std::vector<pid> initialState,
                                   std::vector<pid> finalState,
                                   ComPWA::FourMomentum cmsP4)
    : DalitzKinematics(ParticleStateTransitionKinematicsInfo(
          initialState, finalState, partL, cmsP4, [&finalState]() {
            std::vector<unsigned int> FinalStateEventPositionMapping;
            for (unsigned int i = 0; i < finalState.size(); ++i) {
              FinalStateEventPositionMapping.push_back(i);
            }
            return FinalStateEventPositionMapping;
          }())) {}

DalitzKinematics::DalitzKinematics(
    ParticleStateTransitionKinematicsInfo KinInfo)
    : DalitzKinematics(KinInfo, 1.0) {
  ///  Calculation of n-dimensional phase space volume.
  ///  ToDo: We need to implement an analytical calculation here
}

DalitzKinematics::DalitzKinematics(
    ParticleStateTransitionKinematicsInfo kininfo, double PhspVol)
    : KinematicsInfo(kininfo), PhspVolume(PhspVol),
      M2(kininfo.getInitialStateInvariantMassSquared()) {
  LOG(INFO) << "DalitzKinematics::"
               "DalitzKinematics() | Initialized kinematics "
               "for reaction "
            << KinematicsInfo;
}

double DalitzKinematics::phspVolume() const { return PhspVolume; }

EventCollection
DalitzKinematics::reduceToPhaseSpace(const EventCollection &Events) const {
  EventCollection PhspSample{Events.Pids};

  auto Dataset = convert(Events);
  auto mA = Dataset.Data["mA"];
  auto mB = Dataset.Data["mB"];
  auto mC = Dataset.Data["mC"];
  auto qAB = Dataset.Data["qAB"];
  auto qBC = Dataset.Data["qBC"];
  auto qCA = Dataset.Data["qCA"];

  for (size_t i = 0; i < Events.Events.size(); ++i) {
    double s2 = (qAB[i] + qBC[i] + qCA[i] - mA[i] * mA[i] - mB[i] * mB[i] -
                 mC[i] * mC[i]);

    if (s2 < M2)
      PhspSample.Events.push_back(Events.Events[i]);
  }
  return PhspSample;
} // namespace EvtGen

ComPWA::Data::DataSet
DalitzKinematics::convert(const ComPWA::EventCollection &DataSample) const {

  ComPWA::Data::DataSet Dataset;
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
    throw ComPWA::BadParameter("DalitzKinematics::convert() | number of PIDs "
                               "not equal to number of four-momenta");
  }

  std::vector<double> mA, mB, mC, qAB, qBC, qCA, Weights;
  for (auto const &Event : DataSample.Events) {
    mA.push_back(Event.FourMomenta[0].invariantMass());
    mB.push_back(Event.FourMomenta[1].invariantMass());
    mC.push_back(Event.FourMomenta[2].invariantMass());
    qAB.push_back(
        (Event.FourMomenta[0] + Event.FourMomenta[1]).invariantMass());
    qBC.push_back(
        (Event.FourMomenta[1] + Event.FourMomenta[2]).invariantMass());
    qCA.push_back(
        (Event.FourMomenta[2] + Event.FourMomenta[0]).invariantMass());
    Weights.push_back(Event.Weight);
  }

  Dataset.Data.insert(std::make_pair("mA", mA));
  Dataset.Data.insert(std::make_pair("mB", mB));
  Dataset.Data.insert(std::make_pair("mC", mC));
  Dataset.Data.insert(std::make_pair("qAB", qAB));
  Dataset.Data.insert(std::make_pair("qBC", qBC));
  Dataset.Data.insert(std::make_pair("qCA", qCA));
  Dataset.Data.insert(std::make_pair("weights", Weights));

  return Dataset;
}

} // namespace EvtGen
} // namespace Physics
} // namespace ComPWA
