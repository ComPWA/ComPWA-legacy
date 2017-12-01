// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Containt HelicityDecay class which represents a two-body decay within the
/// helicity formalism.
///
#ifndef HelicityDecay_HPP_
#define HelicityDecay_HPP_

#include <memory>
#include <boost/property_tree/ptree.hpp>

#include "Physics/PartialAmplitude.hpp"
#include "Physics/DecayDynamics/AbstractDynamicalFunction.hpp"
#include "Physics/HelicityFormalism/AmpWignerD.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

///
/// \class HelicityDecay
/// HelicityDecay class represents a two-body decay within the helicity
/// formalism.
///
class HelicityDecay : public ComPWA::Physics::PartialAmplitude {

public:
  HelicityDecay(int dataPos, const SubSystem &sys)
      : DataPosition(dataPos), SubSys(sys){};

  virtual HelicityDecay *clone(std::string newName = "") const {
    auto tmp = new HelicityDecay(*this);
    tmp->setName(newName);
    return tmp;
  }

  /// Factory for HelicityDecay
  static std::shared_ptr<ComPWA::Physics::PartialAmplitude>
  Factory(std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
          const boost::property_tree::ptree &pt);

  static boost::property_tree::ptree Save(std::shared_ptr<PartialAmplitude> obj);

  //======= INTEGRATION/NORMALIZATION ===========

  /// Check of parameters have changed and normalization has to be recalculated
  virtual bool isModified() const;

  virtual double normalization() const;

  //================ EVALUATION =================

  /// Evaluate function without normalization
  std::complex<double> evaluateNoNorm(const DataPoint &point) const {
    std::complex<double> result = coefficient();
    result *= AngularDist->evaluate(point, DataPosition + 1, DataPosition + 2);
    result *= DynamicFcn->evaluate(point, DataPosition);

    assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
    return result;
  };
  //============ SET/GET =================

  virtual void parameters(ParameterList &list);

  virtual void parametersFast(std::vector<double> &list) const {
    PartialAmplitude::parametersFast(list);
    DynamicFcn->GetParametersFast(list);
  }

  /// Update parameters to the values given in \p list
  virtual void updateParameters(const ParameterList &list);

  std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD> wignerD() {
    return AngularDist;
  }

  void setWignerD(
      std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD> w) {
    AngularDist = w;
  }

  std::shared_ptr<ComPWA::Physics::DecayDynamics::AbstractDynamicalFunction>
  dynamicalFunction() {
    return DynamicFcn;
  }

  void setDynamicalFunction(
      std::shared_ptr<ComPWA::Physics::DecayDynamics::AbstractDynamicalFunction>
          f) {
    DynamicFcn = f;
  }

  /// Set position of variables within dataPoint
  virtual void setDataPosition(int pos) { DataPosition = pos; }

  /// Get position of variables within dataPoint
  virtual int dataPosition() const { return DataPosition; }

  virtual void setSubSystem(const SubSystem &sys) { SubSys = sys; }

  virtual SubSystem subSystem() const { return SubSys; }

  //=========== FUNCTIONTREE =================
  virtual bool hasTree() const { return true; }

  virtual std::shared_ptr<FunctionTree>
  tree(std::shared_ptr<Kinematics> kin, const ComPWA::ParameterList &sample,
          const ComPWA::ParameterList &toySample, std::string suffix);

protected:
  /// Position where variables are stored in dataPoint.
  /// We expect to find the invariant mass of the system at @param DataPosition,
  /// cosTheta at @param DataPosition+1 and phi at @param DataPosition+2
  int DataPosition;

  ComPWA::SubSystem SubSys;

  std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD> AngularDist;

  std::shared_ptr<ComPWA::Physics::DecayDynamics::AbstractDynamicalFunction>
      DynamicFcn;
};

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA

#endif
