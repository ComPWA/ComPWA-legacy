// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef ABSTRACTDYNAMICALFUNCTION_HPP
#define ABSTRACTDYNAMICALFUNCTION_HPP

#include <complex>

#include <boost/property_tree/ptree.hpp>

#include "Core/Properties.hpp"
#include "Core/DataPoint.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Spin.hpp"
#include "Core/FunctionTree.hpp"

#include "Physics/DecayDynamics/FormFactor.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayDynamics {

class AbstractDynamicalFunction {

public:
  //============ CONSTRUCTION ==================

  AbstractDynamicalFunction(std::string name = "")
      : Name(name), DaughterMasses(std::pair<double, double>(-999, -999)),
        Current_mass(-999){};

  virtual ~AbstractDynamicalFunction(){};

  //======= INTEGRATION/NORMALIZATION ===========
  /// Check of parameters have changed and normalization has to be recalculatecd
  virtual bool isModified() const = 0;
  
  /// Label as modified/unmodified
  virtual void setModified(bool b) = 0;
  
  //================ EVALUATION =================
  
  virtual std::complex<double> evaluate(const ComPWA::DataPoint &point,
                                        int pos) const = 0;

  //============ SET/GET =================

  virtual void setName(std::string n) { Name = n; }

  virtual std::string name() { return Name; }

  virtual void parameters(ParameterList &list);

  virtual void parametersFast(std::vector<double> &list) const {
    list.push_back(GetMass());
  }

  /// Update parameters to the values given in \p par
  virtual void updateParameters(const ParameterList &par) = 0;

  virtual void SetMassParameter(std::shared_ptr<FitParameter> mass) {
   Mass = mass;
  }

  virtual std::shared_ptr<FitParameter> GetMassParameter() { return Mass; }

  virtual void SetMass(double mass) {Mass->setValue(mass); }

  virtual double GetMass() const { return Mass->value(); }

  virtual void SetDecayMasses(std::pair<double, double> m) {
    DaughterMasses = m;
  }

  virtual std::pair<double, double> GetDecayMasses() const {
    return DaughterMasses;
  }

  virtual void SetDecayNames(std::pair<std::string, std::string> n) {
    DaughterNames = n;
  }

  virtual std::pair<std::string, std::string> GetDecayNames() const {
    return DaughterNames;
  }

  virtual ComPWA::Spin GetSpin() const { return J; }

  virtual void SetSpin(ComPWA::Spin spin) { J = spin; }

  virtual ComPWA::Spin GetOrbitalAngularMomentum() const { return L; }

  virtual void SetOrbitalAngularMomentum(ComPWA::Spin orbitL) { L = orbitL; }

  //=========== FUNCTIONTREE =================

  virtual bool hasTree() const { return false; }

  virtual std::shared_ptr<ComPWA::FunctionTree>
  tree(const ComPWA::ParameterList &sample, int pos,
       std::string suffix = "") = 0;

protected:
  std::string Name;

  /// Masses of daughter particles
  std::pair<double, double> DaughterMasses;

  /// Names of daughter particles
  std::pair<std::string, std::string> DaughterNames;

  /// Resonance mass
  std::shared_ptr<ComPWA::FitParameter> Mass;

  /// Resonance spin
  ComPWA::Spin J;
  /// Orbital Angular Momentum between two daughters in Resonance decay
  ComPWA::Spin L;

  /// Temporary value of mass (used to trigger recalculation of normalization)
  double Current_mass;
};

} // ns::DecayDynamics
} // ns::Physics
} // ns::ComPWA

#endif
