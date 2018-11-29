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
  /// Create DalitzKinematics from inital and final state particle lists.
  /// The lists contain the pid of initial and final state. The position of a
  /// particle in initial or final state list is used later on for
  /// identification.
  DalitzKinematics(
      std::shared_ptr<PartList> partL, const std::vector<pid> &initialState,
      const std::vector<pid> &finalState,
      const ComPWA::FourMomentum &cmsP4 = ComPWA::FourMomentum(0, 0, 0, 0));

  /// Create HelicityKinematics from a boost::property_tree.
  /// The tree is expected to contain something like:
  /// \code
  /// <HelicityKinematics>
  ///  <PhspVolume>1.45</PhspVolume>
  ///  <InitialState>
  ///    <Particle Name='jpsi' Id='0'/>
  ///  </InitialState>
  ///  <FinalState>
  ///    <Particle Name='pi0' Id='1'/>
  ///    <Particle Name='gamma' Id='0'/>
  ///    <Particle Name='pi0' Id='2'/>
  ///  </FinalState>
  /// </HelicityKinematics>
  /// \endcode
  /// The Id is the position of the particle in input data.
  /// \see HelicityKinematics(std::vector<pid> initialState, std::vector<pid>
  /// finalState)
  DalitzKinematics(std::shared_ptr<PartList> partL,
                   const boost::property_tree::ptree &pt);

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
  void convert(const Event &event, DataPoint &point) const;

  /// Check if \p point is within phase space boundaries.
  bool isWithinPhsp(const DataPoint &point) const;

  /// Get number of variables that are added to dataPoint
  virtual size_t numVariables() const { return 6; }

  /// Calculation of helicity angle.
  /// See (Martin and Spearman, Elementary Particle Theory. 1970)
  /// \deprecated Only used as cross-check.
  double helicityAngle(double M, double m, double m2, double mSpec,
                       double invMassSqA, double invMassSqB) const;

  int dataID(const ComPWA::Physics::SubSystem &) {
    return 0;
  }

  virtual unsigned int
  getDataID(const ComPWA::Physics::SubSystem &) const {
    return 0;
  }

protected:
  double _M;

  ///  Calculation of n-dimensional phase space volume.
  ///  ToDo: We need to implement an analytical calculation here
  double calculatePhspVolume() const { return 1.0; }

  /// Add \p newSys to list of SubSystems and return its ID.
  /// In case that this SubSystem is already in the list only the ID is
  /// returned.
  int createIndex();
};

} // namespace EvtGenIF
} // namespace Physics
} // namespace ComPWA

#endif
