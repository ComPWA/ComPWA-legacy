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
/// The basic functionality is the calculatation of the kinematics variables
/// from four-momenta.
/// The variables that are used by a Resonance within a
/// AmpIntensity depend on the SubSystem in which the Resonance appears.
/// For the helicity formalism for each SubSystem the invariant mass \f$m^2\f$,
/// and the helicity angles \f$cos\Theta\f$ and \f$\phi\f$ are calculated
/// (\ref helicityangles).
///
/// The SubSystem is basically defined by the four-momenta of the resonance
/// decay
/// products and the sum of four-momenta of all other final state particles in
/// the decay. Since usually a large number of possible SubSystems can be
/// defined
/// in a particle decay but only few of them are used, we introuce a bookkeeping
/// system (\ref dataID). This increases efficiency (cpu+memory) significantly.
///
/// \section dataID DataID
///
/// lsdafds
/// \see dataID(SubSystem s)
///
/// \section helicityangles Helicity Angles
///
/// The helicity angles \f$cos\Theta\f$ and \f$\phi\f$ are defined as follows:
/// A resonance R originating from another particles decays to two final state
/// particles \f$f_1,f_2\f$.
/// The momentum of the resonance R is shown in the rest frame of the parent
/// decay
/// while the momenta of the final state particles are boosted to the rest frame
/// of the resonance. The helicity angles \f$cos\Theta\f$ and \f$\phi\f$ are
/// then
/// calculated between the momentum direction of \f$f_1\f$ (or \f$f_2\f$) and
/// the
/// momentum direction of the resonance.
/// A (two-dimensional) illustration is given below
/// \image html HelicityAngle.png "Helicity angle"
/// \see convert(const Event &event, dataPoint &point, SubSystem sys,
///                     const std::pair<double, double> limits) const;
///
///
class HelicityKinematics : public ComPWA::Kinematics {
public:
  HelicityKinematics(ParticleStateTransitionKinematicsInfo ki, double phspvol);
  HelicityKinematics(ParticleStateTransitionKinematicsInfo ki);
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
  ///     -# Calculate four-momenta of center-of-mass system, resonance and of
  ///        one final state particle  by summing up the corresponding
  ///        four-momenta given by \p Event.
  ///     -# Boost the four-momentum of the resonance to the center-of-mass
  ///        frame of the parent decay.
  ///     -# Boost the four-momentum of the final state particle in the rest
  ///        frame of the resonance.
  ///     -# Calculate \f$\Theta\f$ and \f$\phi\f$ between those boosted
  ///        momenta.
  std::pair<double, double>
  calculateHelicityAngles(const Event &Event, const SubSystem &SubSys) const;

  /// Calculates the squared invariant mass \f$m^2\f$ of list of final state
  /// particles \p FinalStateIDs. The actual final state four momenta are
  /// extracted from the \p Event.
  double calculateInvariantMassSquared(const Event &Event,
                                       const IndexList &FinalStateIDs) const;

  /// Creates a #DataSet from \p Events.
  /// Calulates all registered kinematic variables for all Events. Kinematic
  /// variables can be registered for example via the #registerSubSystem(const
  /// SubSystem &newSys) method (see also other register methods). In this way
  /// only the variables are calculated that are used by the model.
  ComPWA::Data::DataSet
  convert(const std::vector<ComPWA::Event> &Events) const final;

  /// Returns a subset of \p Events that are within phase space boundaries.
  std::vector<ComPWA::Event>
  reduceToPhaseSpace(const std::vector<ComPWA::Event> &Events) const final;

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

private:
  ParticleStateTransitionKinematicsInfo KinematicsInfo;

  double PhspVolume;

  /// List of subsystems for which invariant mass and angles are calculated
  std::unordered_map<SubSystem, std::pair<std::string, std::string>> Subsystems;

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
