// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Physics/EvtGen/DalitzKinematics.hpp"

#include "Core/Event.hpp"
#include "Core/Logging.hpp"
#include "Core/Particle.hpp"
#include "Data/DataSet.hpp"

#include "ThirdParty/qft++/include/qft++/Vector4.h"

#include <algorithm>
#include <cmath>
#include <numeric>

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
    ParticleStateTransitionKinematicsInfo kininfo)
    : DalitzKinematics(kininfo, 1.0) {
  ///  Calculation of n-dimensional phase space volume.
  ///  ToDo: We need to implement an analytical calculation here
}

DalitzKinematics::DalitzKinematics(
    ParticleStateTransitionKinematicsInfo kininfo, double phspvol)
    : KinematicsInfo(kininfo), PhspVolume(phspvol),
      M2(kininfo.getInitialStateInvariantMassSquared()) {
  LOG(INFO) << "DalitzKinematics::"
               "DalitzKinematics() | Initialized kinematics "
               "for reaction "
            << KinematicsInfo;
}

double DalitzKinematics::phspVolume() const { return PhspVolume; }

std::vector<ComPWA::Event> DalitzKinematics::reduceToPhaseSpace(
    const std::vector<ComPWA::Event> &Events) const {
  std::vector<Event> Evts;

  auto Dataset = convert(Events);
  auto mA = Dataset.Data["mA"];
  auto mB = Dataset.Data["mB"];
  auto mC = Dataset.Data["mC"];
  auto qAB = Dataset.Data["qAB"];
  auto qBC = Dataset.Data["qBC"];
  auto qCA = Dataset.Data["qCA"];

  for (size_t i = 0; i < Events.size(); ++i) {
    double s2 = (qAB[i] + qBC[i] + qCA[i] - mA[i] * mA[i] - mB[i] * mB[i] -
                 mC[i] * mC[i]);

    if (s2 < M2)
      Evts.push_back(Events[i]);
  }
  return Evts;
}

ComPWA::Data::DataSet
DalitzKinematics::convert(const std::vector<Event> &Events) const {

  ComPWA::Data::DataSet Dataset;

  std::vector<double> mA, mB, mC, qAB, qBC, qCA, weights;
  for (auto const &event : Events) {
    mA.push_back(event.FourMomenta[0].invariantMass());
    mB.push_back(event.FourMomenta[1].invariantMass());
    mC.push_back(event.FourMomenta[2].invariantMass());
    qAB.push_back(
        (event.FourMomenta[0] + event.FourMomenta[1]).invariantMass());
    qBC.push_back(
        (event.FourMomenta[1] + event.FourMomenta[2]).invariantMass());
    qCA.push_back(
        (event.FourMomenta[2] + event.FourMomenta[0]).invariantMass());
    weights.push_back(event.Weight);
  }

  Dataset.Data.insert(std::make_pair("mA", mA));
  Dataset.Data.insert(std::make_pair("mB", mB));
  Dataset.Data.insert(std::make_pair("mC", mC));
  Dataset.Data.insert(std::make_pair("qAB", qAB));
  Dataset.Data.insert(std::make_pair("qBC", qBC));
  Dataset.Data.insert(std::make_pair("qCA", qCA));
  Dataset.Data.insert(std::make_pair("weights", weights));

  return Dataset;
}

// int DalitzKinematics::createIndex() {

//   if (VariableNames.size())
//     return 1;

//   VariableNames.push_back("mA");
//   VariableNames.push_back("mB");
//   VariableNames.push_back("mC");
//   VariableNames.push_back("qAB");
//   VariableNames.push_back("qBC");
//   VariableNames.push_back("qCA");
//   //}
//   return 0;
// }

} // namespace EvtGen
} // namespace Physics
} // namespace ComPWA
