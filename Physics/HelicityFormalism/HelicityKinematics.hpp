// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Contains HelicityKinematics class
///

#ifndef PHYSICS_HELICITYFORMALISM_HELICITYKINEMATICS_HPP_
#define PHYSICS_HELICITYFORMALISM_HELICITYKINEMATICS_HPP_

#include <vector>

#include "Core/Kinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

/// \class HelicityKinematics
/// Implementation of the ComPWA::Kinematics interface for amplitude models
/// using the helicity formalism.

class HelicityKinematics : public ComPWA::Kinematics {

public:
  /// Create HelicityKinematics from inital and final state particle lists.
  /// The lists contain the pid of initial and final state. The position of a
  /// particle in initial or final state list is used later on for
  /// identification.
  HelicityKinematics(std::vector<pid> initialState,
                     std::vector<pid> finalState);

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
  HelicityKinematics(boost::property_tree::ptree pt);

  /// Delete copy constructor. For each Kinematics in the analysis only
  /// one instance should exist since Kinematics does the bookkeeping for which
  /// SubSystems variables needs to be calculated. That instance can then be
  /// passed as (smart) pointer.
  /// Note1: Not sure if we also should delete the move constructor.
  /// Note2: We have to delete the copy constructor in Base and Derived classes.
  HelicityKinematics(const HelicityKinematics &that) = delete;

  /// Fill \p point from \p event.
  /// The triple (\f$m^2,cos\Theta, \phi\f$) is added to dataPoint for each
  /// SubSystem that was requested via GetDataID(SubSystem). In this way only
  /// the variables are calculated that are used by the model.
  void EventToDataPoint(const Event &event, dataPoint &point) const;

  /// Fill \p point with variables for \p sys.
  /// The triple (\f$m^2, cos\Theta, \phi\f$) is added to dataPoint for
  /// SubSystem \p sys. Invariant mass limits of the SubSystem can be given via
  /// \p limits.
  void EventToDataPoint(const Event &event, dataPoint &point, SubSystem sys,
                        const std::pair<double, double> limits) const;

  /// Fill \p point with variables for \p sys.
  /// \see EventToDataPoint(const Event &event, dataPoint &point, SubSystem sys,
   ///                     const std::pair<double, double> limits) const;
  void EventToDataPoint(const Event &event, dataPoint &point,
                        SubSystem sys) const;

  /// Check if \p point is within phase space boundaries.
  bool IsWithinPhsp(const dataPoint &point) const;

  virtual bool IsWithinBoxPhsp(int idA, int idB, double varA,
                               double varB) const {
    return true;
  }

  /// Calculation of helicity angle.
  /// Deprecated. Only used as cross-check.
  /// See (Martin and Spearman, Elementary Particle Theory. 1970)
  double HelicityAngle(double M, double m, double m2, double mSpec,
                       double invMassSqA, double invMassSqB) const;

  /// Get ID of data for \p subSys.
  /// In case that the ID was not requested before the subsystem is added to
  /// the list and variables (m^2, cosTheta, phi) are calculated in
  /// #EventToDataPoint()
  /// and added to each dataPoint.
  virtual int GetDataID(const SubSystem subSys) {
    // We calculate the variables currently for two-body decays
    if (subSys.GetFinalStates().size() != 2)
      return 0;
    int pos = createIndex(subSys);
    //    LOG(trace) << " Subsystem " << s << " has dataID " << pos;
    return pos;
  }

  /// Get ID of data for subsystem defined by \p recoilS and \p finalS.
  /// \see GetDataID(SubSystem s)
  virtual int GetDataID(std::vector<int> recoilS, std::vector<int> finalA,
                        std::vector<int> finalB) {
    return GetDataID(SubSystem(recoilS, finalA, finalB));
  }

  /// Get SubSystem from \p pos in list
  virtual SubSystem GetSubSystem(int pos) const {
    return _listSubSystem.at(pos);
  }

  /// Get number of variables that are added to dataPoint
  virtual size_t GetNVars() const { return _listSubSystem.size() * 3; }

  /// Get phase space bounds for the invariant mass of \p subSys.
  virtual const std::pair<double, double> &
  GetInvMassBounds(const SubSystem subSys) const;

  virtual const std::pair<double, double> &GetInvMassBounds(int sysID) const;

protected:
  ///  Calculation of n-dimensional phase space volume.
  ///  ToDo: We need to implement an analytical calculation here
  double calculatePSArea() { return 1.0; }

  /// Invariant mass of decaying particle
  double _M;

  /// Spin of initial state
  ComPWA::Spin _spinM;

  /// List of subsystems for which invariant mass and angles are calculated
  std::vector<SubSystem> _listSubSystem;
  
  /// Invariant mass bounds for each SubSystem
  std::vector<std::pair<double, double>> _invMassBounds;

  std::pair<double, double> CalculateInvMassBounds(const SubSystem sys) const;

  /// Add \p newSys to list of SubSystems and return its ID.
  /// In case that this SubSystem is already in the list only the ID is
  /// returned.
  int createIndex(SubSystem const &newSys) {
    int results =
        std::find(_listSubSystem.begin(), _listSubSystem.end(), newSys) -
        _listSubSystem.begin();
    if (results == _listSubSystem.size()) {
      _listSubSystem.push_back(newSys);
      _invMassBounds.push_back(CalculateInvMassBounds(newSys));
    }
    return results;
  }
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYFORMALISM_HELICITYKINEMATICS_HPP_ */
