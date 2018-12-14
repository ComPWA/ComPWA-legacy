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
  using IndexListTuple = std::tuple<IndexList, IndexList, IndexList, IndexList>;

public:
  HelicityKinematics(const ParticleStateTransitionKinematicsInfo &ki,
                     double phspvol);
  HelicityKinematics(const ParticleStateTransitionKinematicsInfo &ki);
  /// Create HelicityKinematics from inital and final state particle lists.
  /// The lists contain the pid of initial and final state. The position of a
  /// particle in initial or final state list is used later on for
  /// identification.
  HelicityKinematics(
      std::shared_ptr<PartList> partL, const std::vector<pid> &initialState,
      const std::vector<pid> &finalState,
      const ComPWA::FourMomentum &cmsP4 = ComPWA::FourMomentum(0, 0, 0, 0));

  /// Delete copy constructor. For each Kinematics in the analysis only
  /// one instance should exist since Kinematics does the bookkeeping which
  /// SubSystems variables are needs to be calculated. That instance can then be
  /// passed as (smart) pointer.
  /// Note1: Not sure if we also should delete the move constructor.
  /// Note2: We have to delete the copy constructor in Base and Derived classes.
  HelicityKinematics(const HelicityKinematics &that) = delete;

  /// Create \p point from \p event.
  /// For each SubSystem stored via dataID(const SubSystem subSys) function
  /// #convert(const Event&, dataPoint&, SubSystem,
  /// const std::pair<double, double>) is called. In this way only
  /// the variables are calculated that are used by the model.
  DataPoint convert(const Event &event) const;

  /// Fill \p point with variables for \p sys.
  /// The triple (\f$m^2, cos\Theta, \phi\f$) is added to dataPoint for
  /// SubSystem \p sys. Invariant mass limits of the SubSystem can be given via
  /// \p limits.
  /// For a more general introduction see \ref helicityangles. The step-by-step
  /// procedure to calculate the helicity angles is:
  ///     -# Calculate four-momenta of center-of-mass system, resonance and of
  ///        one final state particle  by summing up the corresponding
  ///        four-momenta given by \p event.
  ///     -# The invariant mass \f$m^2\f$ is directly calculated from the
  ///        four-momentum of the resonance.
  ///     -# Boost the four-momentum of the resonance to the center-of-mass
  ///        frame of the parent decay.
  ///     -# Boost the four-momentum of the final state particle in the rest
  ///        frame of the resonance.
  ///     -# Calculate \f$cos\Theta\f$ and \f$\phi\f$ between those boosted
  ///        momenta.
  void convert(const Event &event, DataPoint &point, const SubSystem &sys,
               const std::pair<double, double> limits) const;

  /// Fill \p point with variables for \p sys.
  /// \see convert(const Event &event, dataPoint &point, SubSystem sys,
  ///                     const std::pair<double, double> limits) const;
  void convert(const Event &event, DataPoint &point,
               const SubSystem &sys) const;

  /// Check if \p point is within phase space boundaries.
  bool isWithinPhaseSpace(const DataPoint &point) const final;

  /// Get ID of data for \p subSys.
  unsigned int getDataID(const SubSystem &subSys) const;

  void createAllSubsystems();

  /// Add \p newSys to list of SubSystems and return its ID.
  /// In case that this SubSystem is already in the list only the ID is
  /// returned.
  unsigned int addSubSystem(const SubSystem &newSys);

  /// Add SubSystem from \p pos indices of final state particles
  unsigned int addSubSystem(const std::vector<unsigned int> &FinalA,
                            const std::vector<unsigned int> &FinalB,
                            const std::vector<unsigned int> &Recoil,
                            const std::vector<unsigned int> &ParentRecoil);

  /// Get SubSystem from \p pos in list
  SubSystem subSystem(unsigned int pos) const { return Subsystems.at(pos); }

  /// Get SubSystem from \p pos in list
  std::vector<SubSystem> subSystems() const { return Subsystems; }

  /// Get number of variables that are added to dataPoint
  size_t numVariables() const { return Subsystems.size() * 3; }

  /// Get phase space bounds for the invariant mass of \p subSys.
  const std::pair<double, double> &invMassBounds(const SubSystem &subSys) const;

  const std::pair<double, double> &invMassBounds(int sysID) const;

  /// Calculation of helicity angle.
  /// See (Martin and Spearman, Elementary Particle Theory. 1970)
  /// \deprecated Only used as cross-check.
  double helicityAngle(double M, double m, double m2, double mSpec,
                       double invMassSqA, double invMassSqB) const;

  std::vector<std::string> getKinematicVariableNames() const {
    return VariableNames;
  }

  double phspVolume() const;

  const ParticleStateTransitionKinematicsInfo &
  getParticleStateTransitionKinematicsInfo() const {
    return KinematicsInfo;
  }

private:
  ParticleStateTransitionKinematicsInfo KinematicsInfo;

  double PhspVolume;

  /// List of subsystems for which invariant mass and angles are calculated
  std::vector<SubSystem> Subsystems;

  /// Invariant mass bounds for each SubSystem
  std::vector<std::pair<double, double>> InvMassBounds;

  std::vector<std::string> VariableNames;

  std::pair<double, double> calculateInvMassBounds(const SubSystem &sys) const;

  IndexListTuple sortSubsystem(const IndexListTuple &SubSys) const;
  std::vector<std::pair<IndexList, IndexList>>
  redistributeIndexLists(const IndexList &A, const IndexList &B) const;
};

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA

#endif
