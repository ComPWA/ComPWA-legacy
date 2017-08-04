//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

#ifndef PHYSICS_HELICITYAMPLITUDE_HELICITYKINEMATICS_HPP_
#define PHYSICS_HELICITYAMPLITUDE_HELICITYKINEMATICS_HPP_

#include <vector>

#include "Core/Kinematics.hpp"
#include "Physics/HelicityFormalism/SubSystem.hpp"

//class Particle;

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {


/*! HelicityKinematics class.
   Implementation of the ComPWA::Kinematics interface for amplitude models
   using the helicity formalism.
 */
class HelicityKinematics : public ComPWA::Kinematics {

public:

  /*! Create HelicityKinematics from inital and final state particle lists.
   * The lists contain the pid of initial and final state. The position of a
   * particle in initial or final state list is used later on for
   * identification.
   */
  HelicityKinematics(std::vector<pid> initialState,
                     std::vector<pid> finalState);
  
  /*! Create HelicityKinematics from a property_tree.
   * The tree is expected to contain something like:
   * @code
    <HelicityKinematics>
      <PhspVolume>1.45</PhspVolume>
      <InitialState>
        <Particle Name='jpsi' Id='0'/>
      </InitialState>
      <FinalState>
        <Particle Name='pi0' Id='1'/>
        <Particle Name='gamma' Id='0'/>
        <Particle Name='pi0' Id='2'/>
      </FinalState>
    </HelicityKinematics>
    @endcode
   * The Id is the position of the particle in input data.
   * @see HelicityKinematics(std::vector<pid> initialState, std::vector<pid>
   finalState)
   */
  HelicityKinematics(boost::property_tree::ptree pt);

  /*! Fill dataPoint point.
   * The triple (\f$m^2,cos\Theta, \phi\f$) is added to dataPoint for each
   * SubSystem that was requested via GetDataID(SubSystem).
   */
  void EventToDataPoint(const Event &event, dataPoint &point) const;

  /*! Fill dataPoint point with variables for SubSystem.
   * The triple (\f$m^2,cos\Theta, \phi\f$) is added to dataPoint for
   * SubSystem sys.
   */
  void EventToDataPoint(const Event &event, dataPoint &point, SubSystem sys,
                        const std::pair<double, double> limits) const;

  void EventToDataPoint(const Event &event, dataPoint &point,
                        SubSystem sys) const;
  
  /*! Check if dataPoint is within phase space boundaries.
   */
  bool IsWithinPhsp(const dataPoint &point) const;

  virtual bool IsWithinBoxPhsp(int idA, int idB, double varA,
                               double varB) const {
    return true;
  }

  //! get mass of particles
  virtual double GetMass(unsigned int num) const { return 0; }

  //! get mass of particles
  virtual double GetMass(std::string name) const { return 0; }

  //! Get mass of mother particle
  virtual double GetMotherMass() const { return _M; }

  //! Get number of particles
  virtual std::size_t GetNumberOfParticles() const {
    return _finalState.size();
  }
  
  double HelicityAngle(double M, double m, double m2, double mSpec,
                       double invMassSqA, double invMassSqB) const;

  /*! Get ID of data for SubSystem #s.
   * Incase that the ID was not requested before the subsystem is added to the
   * list and variables (m^2, cosTheta, phi) are calculated in
   * #EventToDataPoint()
   * and added to each dataPoint.
   */
  virtual int GetDataID(const SubSystem s) {
    //We calculate the variables currently for two-body decays
    if( s.GetFinalStates().size() != 2) return 0;
    int pos = createIndex(s);
    //    LOG(trace) << " Subsystem " << s << " has dataID " << pos;
    return pos;
  }

  /*! Get ID of data for subsystem defined by recoilS and finalS.
   * @see GetDataID(SubSystem s)
   */
  virtual int GetDataID(std::vector<int> recoilS, std::vector<int> finalA,
                        std::vector<int> finalB) {
    return GetDataID(SubSystem(recoilS, finalA, finalB));
  }

  //! Get SubSystem from pos in list
  virtual SubSystem GetSubSystem(int pos) const {
    return _listSubSystem.at(pos);
  }

  //! Get number of variables that are added to dataPoint
  virtual size_t GetNVars() const { return _listSubSystem.size() * 3; }

  /*! Get phase space bounds for the invariant mass of SubSystem sys.
   */
  virtual const std::pair<double, double> &
  GetInvMassBounds(const SubSystem sys) const;
  
  virtual const std::pair<double, double> &GetInvMassBounds(int sysID) const;

protected:
  /*! Calculation of n-dimensional phase space volume.
   * We calculate the phase space volume using an estimate of the (constant)
   * event density in phase space. We follow this procedure:
   *   1. Generate phase space sample
   *   2. For each event calculate the distance to each other event using
   *      #EventDistance and store it in a vector
   *   3. Sort the vector and select the distance within which
   *      #localDensityNumber of events are found
   *   4. Calculate the volume of a n-dimensional sphere with radius according
   *      to 3.
   *   5. Calculate local density using #localDensityNumber and the volume
   *   6. Calculate the average density. At this point we assume that points
   *      are distributed uniformly across phase space. Therefore, the set
   *      #_irreducibleSetOfVariables needs to represent a flat phase space.
   *   7. From the event density in phase space and the number of events in
   * the
   *      sample the volume can be calculated.
   */
  double calculatePSArea();

  /*! Calculate distance in phase space of two events.
   * We calculate the eucleadean distance of the variables listed in
   * #_irreducibleSetOfVariables for both events. This procedure does not make
   * sense if the set is not irreducible.
   */
  double EventDistance(Event &evA, Event &evB) const;

  //! Invariant mass
  double _M;

  //! Invariant mass squared
  double _Msq;

  //! Spin of initial state
  ComPWA::Spin _spinM;

  //! List of subsystems for which invariant mass and angles are calculated
  std::vector<SubSystem> _listSubSystem;
  std::vector<std::pair<double, double>> _invMassBounds;

  std::pair<double, double> CalculateInvMassBounds(const SubSystem sys) const;

  int createIndex(SubSystem const &newValue) {
    int results =
        std::find(_listSubSystem.begin(), _listSubSystem.end(), newValue) -
        _listSubSystem.begin();
    if (results == _listSubSystem.size()) {
      _listSubSystem.push_back(newValue);
      _invMassBounds.push_back(CalculateInvMassBounds(newValue));
    }
    return results;
  }
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_HELICITYKINEMATICS_HPP_ */
