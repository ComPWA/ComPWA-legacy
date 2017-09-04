// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Containt PartialDecay class which represents a two-body decay within the
/// helicity formalism.
///
#ifndef PARTIALDECAY_HPP_
#define PARTIALDECAY_HPP_

#include <memory>
#include <boost/property_tree/ptree.hpp>

#include "Physics/Resonance.hpp"
#include "Physics/DecayDynamics/AbstractDynamicalFunction.hpp"
#include "Physics/HelicityFormalism/AmpWignerD.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

///
/// \class PartialDecay
/// PartialDecay class represents a two-body decay within the helicity
/// formalism.
///
class PartialDecay : public ComPWA::Physics::Resonance {

public:
  PartialDecay(SubSystem sys) : _subSystem(sys) {};
  
  virtual PartialDecay *Clone(std::string newName = "") const {
    auto tmp = new PartialDecay(*this);
    tmp->SetName(newName);
    return tmp;
  }

  /// Factory for PartialDecay
  static std::shared_ptr<ComPWA::Physics::Resonance>
  Factory(std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
          const boost::property_tree::ptree &pt);

  static boost::property_tree::ptree Save(std::shared_ptr<Resonance> obj);

  //======= INTEGRATION/NORMALIZATION ===========

  /// Check of parameters have changed and normalization has to be recalculated
  virtual bool CheckModified() const;

  virtual double GetNormalization() const;

  //================ EVALUATION =================

  /// Evaluate function without normalization
  std::complex<double> EvaluateNoNorm(const dataPoint &point) const {
    std::complex<double> result = GetCoefficient();
    result *= _angD->Evaluate(point, _dataPos + 1, _dataPos + 2);
    result *= _dynamic->Evaluate(point, _dataPos);

    assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
    return result;
  };
  //============ SET/GET =================

  virtual void GetParameters(ParameterList &list);

  virtual void GetParametersFast(std::vector<double> &list) const {
    Resonance::GetParametersFast(list);
    _dynamic->GetParametersFast(list);
  }
  
  /// Update parameters to the values given in \p list
  virtual void UpdateParameters(const ParameterList &list);

  std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD> GetWignerD() {
    return _angD;
  }

  void SetWignerD(
      std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD> w) {
    _angD = w;
  }

  std::shared_ptr<ComPWA::Physics::DecayDynamics::AbstractDynamicalFunction>
  GetDynamicalFunction() {
    return _dynamic;
  }

  void SetDynamicalFunction(
      std::shared_ptr<ComPWA::Physics::DecayDynamics::AbstractDynamicalFunction>
          f) {
    _dynamic = f;
  }

  /// Set position of variables within dataPoint
  virtual void SetDataPosition(int pos) { _dataPos = pos; }

  /// Get position of variables within dataPoint
  virtual int GetDataPosition() const { return _dataPos; }

  virtual void SetSubSystem(SubSystem sys) { _subSystem = sys; }

  virtual SubSystem GetSubSystem() const { return _subSystem; }

  //=========== FUNCTIONTREE =================
  virtual bool HasTree() const { return true; }

  virtual std::shared_ptr<FunctionTree>
  GetTree(std::shared_ptr<Kinematics> kin, const ComPWA::ParameterList &sample,
          const ComPWA::ParameterList &toySample, std::string suffix);

protected:
  /// Position where variables are stored in dataPoint.
  /// We expect to find the invariant mass of the system at @param _dataPos,
  /// cosTheta at @param _dataPos+1 and phi at @param _dataPos+2
  int _dataPos;

  ComPWA::SubSystem _subSystem;

  std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD> _angD;

  std::shared_ptr<ComPWA::Physics::DecayDynamics::AbstractDynamicalFunction>
      _dynamic;
};

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA

#endif
