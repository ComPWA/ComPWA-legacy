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
#include "Physics/qft++/Vector4.h"

class Particle;

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

static const char *formFactorTypeString[] = {"noFormFactor", "BlattWeisskopf",
                                             "CrystalBarrel"};

enum formFactorType { noFormFactor = 0, BlattWeisskopf = 1, CrystalBarrel = 2 };

template <typename T>
int createIndex(std::vector<T> &references, T const &newValue) {
  int results = std::find(references.begin(), references.end(), newValue) -
                references.begin();
  if (results == references.size()) {
    references.push_back(newValue);
  }
  return results;
}

/*! Definition of a two-body decay node within a sequential decay tree.
 Class contains lists for both final states of the two-body decay and a list 
 for all recoiling particles. This information is needed to calculate
 invariant mass and angles at a two-body decay node.
 */
class SubSystem {
public:
  SubSystem(){};

  SubSystem(std::vector<int> recoilS, std::vector<int> finalA,
            std::vector<int> finalB)
      : _recoilState(recoilS), _finalStateA(finalA), _finalStateB(finalB) {

    // Creating unique title
    std::stringstream stream;
    stream << "(";
    for (auto i : _recoilState)
      stream << std::to_string(i);
    stream << ")->(";
    for (auto i : _finalStateA)
      stream << std::to_string(i);
    stream << ")+(";
    for (auto i : _finalStateB)
      stream << std::to_string(i);
    stream << ")";
    title = stream.str();
    // LOG(trace) << "SubSystem::SubSystem() | Creating sub system "<<title;
  }

  bool operator==(const SubSystem &b) const {
    if (_recoilState == b._recoilState && _finalStateA == b._finalStateA &&
        _finalStateB == b._finalStateB)
      return true;
    return false;
  }

  friend std::ostream &operator<<(std::ostream &stream, const SubSystem &s) {
    stream << "( ";
    for (auto i : s._recoilState)
      stream << i << " ";
    stream << ") <-> ( ";
    for (auto i : s._finalStateA)
      stream << i << " ";
    stream << ")( ";
    for (auto i : s._finalStateB)
      stream << i << " ";
    stream << ")";
    return stream;
  }

  std::vector<int> GetRecoilState() const { return _recoilState; }
  std::vector<int> GetFinalStateA() const { return _finalStateA; }
  std::vector<int> GetFinalStateB() const { return _finalStateB; }
  void SetRecoilState(std::vector<int> r) { _recoilState = r; }
  void SetFinalStateA(std::vector<int> f) { _finalStateA = f; }
  void SetFinalStateB(std::vector<int> f) { _finalStateB = f; }

protected:
  std::string title;
  std::vector<int> _recoilState;
  std::vector<int> _finalStateA;
  std::vector<int> _finalStateB;
};

/*! HelicityKinematics class.
   Implementation of the ComPWA::Kinematics interface for amplitude models
   using the helicity formalism.
 */
class HelicityKinematics : public ComPWA::Kinematics {

public:
  /*! Create instance of HelicityKinematics.
   * @see HelicityKinematics(std::vector<pid>, std::vector<pid>);
   */
  static Kinematics *CreateInstance(std::vector<int> initialState,
                                    std::vector<int> finalStateIds) {
    if (_inst) {
      throw std::runtime_error(
          "HelicityKinematics::createInstance() | Instance already exists!");
    }
    _inst = new HelicityKinematics(initialState, finalStateIds);
    return _inst;
  }

  /*! Create instance of HelicityKinematics.
   * @see HelicityKinematics(std::vector<pid>, std::vector<pid>);
   */
  static Kinematics *CreateInstance(boost::property_tree::ptree pt) {
    if (_inst) {
      throw std::runtime_error(
          "HelicityKinematics::createInstance() | Instance already exists!");
    }
    _inst = new HelicityKinematics(pt);
    return _inst;
  }

  /*! Fill dataPoint point.
   * The triple (\f$m^2,cos\Theta, \phi\f$) is added to dataPoint for each
   * SubSystem that was requested via GetDataID(SubSystem).
   */
  void EventToDataPoint(const Event &event, dataPoint &point) const;

  // delete methods to ensure that there will only be one instance
  HelicityKinematics() = delete;

  HelicityKinematics(const HelicityKinematics &) = delete;

  void operator=(const HelicityKinematics &) = delete;

  //! Check if dataPoint is within phase space boundaries.
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
  virtual unsigned int GetNumberOfParticles() const {
    return _finalState.size();
  }

  /*! Get ID of data for SubSystem #s.
   * Incase that the ID was not requested before the subsystem is added to the
   * list and variables (m^2, cosTheta, phi) are calculated in
   * #EventToDataPoint()
   * and added to each dataPoint.
   */
  virtual int GetDataID(SubSystem s) {
    int pos = createIndex(_listSubSystem, s);
    LOG(trace) << " Subsystem " << s << " has dataID " << pos;
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
  virtual unsigned int GetNVars() const { return _listSubSystem.size() * 3; }

  /*! Get phase space bounds for the invariant mass of SubSystem sys.
   */
  virtual std::pair<double, double> GetInvMassBounds(SubSystem sys) const;

  //! Calculate form factor
  static double
  FormFactor(double sqrtS, double ma, double mb, double spin,
             double mesonRadius,
             formFactorType type = formFactorType::BlattWeisskopf);

  //! Calculate form factor
  static double
  FormFactor(double sqrtS, double ma, double mb, double spin,
             double mesonRadius, std::complex<double> qValue,
             formFactorType type = formFactorType::BlattWeisskopf);

  /*! Specify an irreducible set of variables.
   * An irreducible set of variables that represent a phase space position.
   * @see _irreducibleSetOfVariables
   */
  virtual void SetIrreducibleSetOfVariables(std::vector<int> set) {
    _irreducibleSetOfVariables = set;
  }

  /*! Get an irreducible set of variables.
   * @see SetIrreducibleSetOfVariables(std::vector<int> set)
   * @see _irreducibleSetOfVariables
   */
  virtual std::vector<int> GetIrreducibleSetOfVariables() const {
    return _irreducibleSetOfVariables;
  }

protected:
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
      <InitialState>
        <Particle Name='jpsi' />
      </InitialState>
      <FinalState>
        <Particle Name='pi0' />
        <Particle Name='gamma' />
        <Particle Name='pi0' />
      </FinalState>
    </HelicityKinematics>
    @endcode
   * @see HelicityKinematics(std::vector<pid> initialState, std::vector<pid>
   finalState)
   */
  HelicityKinematics(boost::property_tree::ptree pt);

  virtual ~HelicityKinematics();

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
   *   7. From the event density in phase space and the number of events in the
   *      sample the volume can be calculated.
   */
  double calculatePSArea();

  /*! Calculate distance in phase space of two events.
   * We calculate the eucleadean distance of the variables listed in
   * #_irreducibleSetOfVariables for both events. This procedure does not make
   * sense if the set is not irreducible.
   */
  double EventDistance(Event &evA, Event &evB) const;

  /*! Irreducible set of variable that represent a position in phase space.
   * Listed are the positions in dataPoint. The size of the set corresponds to
   * the degrees of freedom of the decay.
   * It is important that phase space is distributed uniformly in those
   * variables otherwise the normalization of intensities is not correct.
   */
  std::vector<int> _irreducibleSetOfVariables;

  //! Invariant mass
  double _M;

  //! Invariant mass squared
  double _Msq;

  //! Spin of initial state
  ComPWA::Spin _spinM;

  //! List of subsystems for which invariant mass and angles are calculated
  std::vector<SubSystem> _listSubSystem;
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_HELICITYKINEMATICS_HPP_ */
