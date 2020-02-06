// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Contains DalitzKinematics class.
///

#ifndef DALITZKINEMATICS_HPP_
#define DALITZKINEMATICS_HPP_

#include <vector>

#include "Core/Kinematics.hpp"
#include "Physics/ParticleStateTransitionKinematicsInfo.hpp"
#include "Physics/SubSystem.hpp"

namespace ComPWA {
namespace Physics {
namespace EvtGen {

///
/// \class DalitzKinematics
/// Implementation of the ComPWA::Kinematics interface for amplitude models
/// using the helicity formalism.
/// The basic functionality is the calculatation of the kinematics variables
/// from four-momenta.

class DalitzKinematics : public ComPWA::Kinematics {

public:
  DalitzKinematics(ParticleStateTransitionKinematicsInfo kininfo,
                   double phspvol);
  DalitzKinematics(ParticleStateTransitionKinematicsInfo kininfo);

  /// Create DalitzKinematics from inital and final state particle lists.
  /// The lists contain the pid of initial and final state. The position of a
  /// particle in initial or final state list is used later on for
  /// identification.
  DalitzKinematics(ComPWA::ParticleList partL, std::vector<pid> initialState,
                   std::vector<pid> finalState,
                   ComPWA::FourMomentum cmsP4 = ComPWA::FourMomentum(0, 0, 0,
                                                                     0));

  /// Delete copy constructor. For each Kinematics in the analysis only
  /// one instance should exist since Kinematics does the bookkeeping which
  /// SubSystems variables are needs to be calculated. That instance can then be
  /// passed as (smart) pointer.
  DalitzKinematics(const DalitzKinematics &that) = delete;
  DalitzKinematics(DalitzKinematics &&that) = default;

  ComPWA::Data::DataSet convert(const std::vector<Event> &Events) const final;

  /// Returns a subset of \p Events that are within phase space boundaries.
  std::vector<ComPWA::Event>
  reduceToPhaseSpace(const std::vector<ComPWA::Event> &Events) const final;

  double phspVolume() const;

private:
  ParticleStateTransitionKinematicsInfo KinematicsInfo;

  double PhspVolume;

  double M2;
};

} // namespace EvtGen
} // namespace Physics
} // namespace ComPWA

#endif
