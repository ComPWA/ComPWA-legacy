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

//#include "Physics/HelicityAmplitude/ParticleStateDefinitions.hpp"
//#include "Physics/HelicityAmplitude/FinalStateParticleCombinatorics.hpp"

#include "Physics/qft++/Vector4.h"

class Particle;

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

template <typename T>
int createIndex(std::vector<T> &references, T const &newValue) {
  int results = std::find(references.begin(), references.end(), newValue) -
                references.begin();
  if (results == references.size()) {
    references.push_back(newValue);
  }
  return results;
}

class SubSystem {
public:
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
    stream << "recoilParticles -> (";
    for (auto i : s._recoilState)
      stream << i << " ";
    stream << ") finalParticlesA -> (";
    for (auto i : s._finalStateA)
      stream << i << " ";
    stream << ") finalParticlesB -> (";
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

class HelicityKinematics : public ComPWA::Kinematics {

public:
  static Kinematics *CreateInstance(std::vector<int> initialState,
                                    std::vector<int> finalStateIds) {
    if (_inst) {
      throw std::runtime_error(
          "HelicityKinematics::createInstance() | Instance already exists!");
    }
    _inst = new HelicityKinematics(initialState, finalStateIds);
    return _inst;
  }
  
  static Kinematics *CreateInstance(boost::property_tree::ptree pt) {
    if (_inst) {
      throw std::runtime_error(
          "HelicityKinematics::createInstance() | Instance already exists!");
    }
    _inst = new HelicityKinematics(pt);
    return _inst;
  }

  //! converts Event to dataPoint
  void EventToDataPoint(const Event &event, dataPoint &point) const;

  // delete methods to ensure that there will only be one instance
  HelicityKinematics() = delete;

  HelicityKinematics(const HelicityKinematics &) = delete;

  void operator=(const HelicityKinematics &) = delete;

  bool IsWithinPhsp(const dataPoint &point) const;

  virtual bool IsWithinBoxPhsp(int idA, int idB, double varA,
                               double varB) const {
    return true;
  }

  //! get mass of particles
  virtual double GetMass(unsigned int num) const { return 0; }

  //! get mass of particles
  virtual double GetMass(std::string name) const { return 0; }

  //! Get name of mother particle
  virtual std::string GetMotherName() const { return _nameMother; };

  //! Get mass of mother particle
  virtual double GetMotherMass() const { return _M; }

  //! Get number of particles
  virtual unsigned int GetNumberOfParticles() const {
    return _finalState.size();
  }

  //! Get ID of data for SubSystem s
  virtual int GetDataID(SubSystem s) {
    int pos = createIndex(_listSubSystem, s);
    LOG(trace) << " Subsystem " << s << " has dataID " << pos;
    return pos;
  }

  //! Get ID of data for subsystem defined by recoilS and finalS
  virtual int GetDataID(std::vector<int> recoilS, std::vector<int> finalA,
                        std::vector<int> finalB) {
    return GetDataID(SubSystem(recoilS, finalA, finalB));
  }

  //! Get ID of data for SubSystem s
  virtual SubSystem GetSubSystem(int pos) const {
    return _listSubSystem.at(pos);
  }

  //! Get number of variables
  virtual unsigned int GetNVars() const { return _listSubSystem.size() * 3; }

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

protected:
  HelicityKinematics(std::vector<int> initialState,
                     std::vector<int> finalState);

  HelicityKinematics(boost::property_tree::ptree pt);
  
  virtual ~HelicityKinematics();

  double calculatePSArea();

  // Parameters of decaying mother particle (we assume that we have a decay)
  std::string _nameMother; //! name of mother particle
  int _idMother;           //! name of mother particle
  double _M;               //! mass of mother particle
  double _Msq;             //! mass of mother particle
  ComPWA::Spin _spinM;     //! spin of mother particle

  // list of subsystems for which invariant mass and angles are calculated
  std::vector<SubSystem> _listSubSystem;
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_HELICITYKINEMATICS_HPP_ */
