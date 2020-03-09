// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PHYSICS_HELICITYFORMALISM_HELICITYKINEMATICS_HPP_
#define PHYSICS_HELICITYFORMALISM_HELICITYKINEMATICS_HPP_

#include <vector>

#include "Core/Kinematics.hpp"
#include "Physics/ParticleStateTransitionKinematicsInfo.hpp"
#include "Physics/SubSystem.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

///
/// \class HelicityKinematics
/// Implementation of the ComPWA::Kinematics interface for amplitude models
/// using the helicity formalism.
/// The basic functionality is the calculation of the kinematics variables
/// from four-momenta.
/// \see ComPWA::Data::DataSet convert(const EventCollection &Events)
/// const;
///
/// Each SubSystem defines three kinematic variables: the invariant mass
/// \f$m^2\f$, and the helicity angles \f$\Theta\f$ and \f$\phi\f$ are
/// calculated
///
/// A SubSystem uniquely defines a two body decay based on the participating
/// four-momenta: the two final states (which make up the decaying state), the
/// decaying state recoil system, and the parents recoil).
/// Since usually a large number of possible SubSystems can be defined in a
/// particle decay but not all of them are used, bookkeeping system is
/// introduced. This increases efficiency (cpu+memory) significantly.
/// Kinematic variables can be registered for calculation via the register
/// methods:
/// -# \see registerSubSystem();
/// -# \see std::string registerInvariantMassSquared(IndexList System);
/// -# \see std::pair<std::string, std::string> registerHelicityAngles(SubSystem
/// System);
///
class HelicityKinematics : public ComPWA::Kinematics {
public:
  HelicityKinematics(ParticleStateTransitionKinematicsInfo KinInfo,
                     double PhspVol);
  HelicityKinematics(ParticleStateTransitionKinematicsInfo KinInfo);
  /// Create HelicityKinematics from inital and final state particle lists.
  /// The lists contain the pid of initial and final state. The position of a
  /// particle in initial or final state list is used later on for
  /// identification.
  HelicityKinematics(ComPWA::ParticleList partL, std::vector<pid> initialState,
                     std::vector<pid> finalState,
                     ComPWA::FourMomentum cmsP4 = ComPWA::FourMomentum(0, 0, 0,
                                                                       0));

  /// Delete copy constructor. For each Kinematics in the analysis only
  /// one instance should exist since Kinematics does the bookkeeping which
  /// SubSystems variables are needs to be calculated. That instance can then be
  /// passed by reference.
  HelicityKinematics(const HelicityKinematics &that) = delete;
  // explicitly default the move constructor because copy constructor was
  // explicitly deleted and that prevents the automatic creation
  HelicityKinematics(HelicityKinematics &&that) = default;

  /// Calculates the pair of values \f$(\Theta, \phi)\f$ of the Event \p Event
  /// for SubSystem \p SubSys. The step-by-step procedure to calculate the
  /// helicity angles is:
  ///     -# Calculate the four-momentum of center-of-mass system of the decay,
  ///     by summing up all final state particle four-momenta.
  ///     -# Boost the the final state particles, the recoil and the parent
  ///     recoil into the CMS of the decaying state.
  ///     -# Rotate the whole CMS system so that the recoil of the decaying
  ///     state points in the -z-axis direction. This makes the z-axis the path
  ///     of flight of the decaying state.
  ///     -# Then rotate the CMS system so that the parent recoil lies in the
  ///     x-z plane. This defines the angle \f$\phi\f$ of the current
  ///     transformed final state particles as the angle difference with respect
  ///     to the production plane.
  ///     -# The helicity angles \f$\Theta\f$ and \f$\phi\f$ can now simply be
  ///     read of the momenta of the final state particles.
  ///
  /// A (two-dimensional) illustration is given below \image html
  /// HelicityAngle.png "Helicity angle"
  std::pair<double, double>
  calculateHelicityAngles(const Event &Event, const SubSystem &SubSys) const;

  /// Calculates the squared invariant mass \f$m^2\f$ of list of final state
  /// particles \p FinalStateIDs. The actual final state four momenta are
  /// extracted from the \p Event.
  double calculateInvariantMassSquared(const Event &Event,
                                       const IndexList &FinalStateIDs) const;

  /// Creates a DataSet from \p Events.
  /// Calculates all registered kinematic variables for all Events. Kinematic
  /// variables can be registered for example via the #registerSubSystem(const
  /// SubSystem &newSys) method (see also other register methods). In this way
  /// only the variables are calculated that are used by the model.
  ComPWA::Data::DataSet convert(const EventCollection &Events) const final;

  /// Returns a subset of \p Events that are within phase space boundaries.
  EventCollection reduceToPhaseSpace(const EventCollection &Events) const final;

  std::string registerInvariantMassSquared(IndexList System);
  std::pair<std::string, std::string> registerHelicityAngles(SubSystem System);

  void createAllSubsystems();

  /// Add \p NewSys to list of SubSystems and return a tuple of names, that id
  /// the registered kinematic variables. In case that this SubSystem is already
  /// in the list only the variable names is returned.
  std::tuple<std::string, std::string, std::string>
  registerSubSystem(const SubSystem &NewSys);

  /// Add SubSystem from \p pos indices of final state particles
  std::tuple<std::string, std::string, std::string>
  registerSubSystem(const std::vector<unsigned int> &FinalA,
                    const std::vector<unsigned int> &FinalB,
                    const std::vector<unsigned int> &Recoil,
                    const std::vector<unsigned int> &ParentRecoil);

  /// Get phase space bounds for the registered invariant mass with name
  /// \p InvariantMassName.
  const std::pair<double, double> &
  getInvariantMassBounds(const std::string &InvariantMassName) const;

  double phspVolume() const;

  const ParticleStateTransitionKinematicsInfo &
  getParticleStateTransitionKinematicsInfo() const {
    return KinematicsInfo;
  }

  const std::vector<pid> &getFinalStatePIDs() const override {
    return KinematicsInfo.getFinalStatePIDs();
  }

private:
  ParticleStateTransitionKinematicsInfo KinematicsInfo;

  double PhspVolume;

  /// Mapping of subsystems to the corresponding helicity angle variable names
  /// (theta, phi)
  std::unordered_map<SubSystem, std::pair<std::string, std::string>> Subsystems;

  /// Mapping of final state particle index lists to invariant mass variable
  /// name
  std::unordered_map<IndexList, std::string> InvariantMassesSquared;

  /// Invariant mass bounds for each SubSystem
  std::unordered_map<std::string, std::pair<double, double>> InvMassBounds;

  std::pair<double, double>
  calculateInvMassBounds(const IndexList &FinalStateIDs) const;
};

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA

#endif
