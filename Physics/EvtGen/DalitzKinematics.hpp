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
  DalitzKinematics(const ParticleStateTransitionKinematicsInfo &kininfo,
                   double phspvol);
  DalitzKinematics(const ParticleStateTransitionKinematicsInfo &kininfo);

  /// Create DalitzKinematics from inital and final state particle lists.
  /// The lists contain the pid of initial and final state. The position of a
  /// particle in initial or final state list is used later on for
  /// identification.
  DalitzKinematics(
      std::shared_ptr<PartList> partL, const std::vector<pid> &initialState,
      const std::vector<pid> &finalState,
      const ComPWA::FourMomentum &cmsP4 = ComPWA::FourMomentum(0, 0, 0, 0));

  /// Delete copy constructor. For each Kinematics in the analysis only
  /// one instance should exist since Kinematics does the bookkeeping which
  /// SubSystems variables are needs to be calculated. That instance can then be
  /// passed as (smart) pointer.
  /// Note1: Not sure if we also should delete the move constructor.
  /// Note2: We have to delete the copy constructor in Base and Derived classes.
  DalitzKinematics(const DalitzKinematics &that) = delete;

  /// Fill \p point from \p event.
  /// For each SubSystem stored via dataID(const SubSystem subSys) function
  /// #convert(const Event&, dataPoint&, SubSystem,
  /// const std::pair<double, double>) is called. In this way only
  /// the variables are calculated that are used by the model.
  DataPoint convert(const Event &event) const;

  /// Check if \p point is within phase space boundaries.
  bool isWithinPhsp(const DataPoint &point) const;

  /// Get number of variables that are added to dataPoint
  virtual size_t numVariables() const { return 6; }

  /// Calculation of helicity angle.
  /// See (Martin and Spearman, Elementary Particle Theory. 1970)
  /// \deprecated Only used as cross-check.
  double helicityAngle(double M, double m, double m2, double mSpec,
                       double invMassSqA, double invMassSqB) const;

  int dataID(const ComPWA::Physics::SubSystem &) { return 0; }

  virtual unsigned int getDataID(const ComPWA::Physics::SubSystem &) const {
    return 0;
  }

  std::vector<std::string> getKinematicVariableNames() const {
    return VariableNames;
  }

  double phspVolume() const;

private:
  ParticleStateTransitionKinematicsInfo KinematicsInfo;

  double PhspVolume;

  double _M;

  std::vector<std::string> VariableNames;

  /// Add \p newSys to list of SubSystems and return its ID.
  /// In case that this SubSystem is already in the list only the ID is
  /// returned.
  int createIndex();
};

} // namespace EvtGen
} // namespace Physics
} // namespace ComPWA

#endif
