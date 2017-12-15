// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PHYSICS_DECAYDYNAMICS_RELATIVISTICBREITWIGNER_HPP_
#define PHYSICS_DECAYDYNAMICS_RELATIVISTICBREITWIGNER_HPP_

#include <vector>
#include <memory>
#include <boost/property_tree/ptree.hpp>

#include "Core/Spin.hpp"
#include "Core/Functions.hpp"
#include "Core/Exceptions.hpp"
#include "Physics/DecayDynamics/AbstractDynamicalFunction.hpp"
#include "Physics/HelicityFormalism/AmpWignerD.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayDynamics {

class HelicityDecay;

/// \class RelativisticBreitWigner
/// Relativistic Breit-Wigner model with barrier factors
/// The dynamical function implemented here is taken from PDG2014 (Eq.47-22) for
/// the one channel case.
class RelativisticBreitWigner
    : public ComPWA::Physics::DecayDynamics::AbstractDynamicalFunction {

public:
  RelativisticBreitWigner(std::string name, std::pair<std::string,std::string> daughters,
               std::shared_ptr<ComPWA::PartList> partL);

  /// Check of parameters have changed and normalization has to be recalculatecd
  virtual bool isModified() const;
   
  //================ EVALUATION =================
  std::complex<double> evaluate(const ComPWA::DataPoint &point, int pos) const;

  /// Dynamical Breit-Wigner function.
  /// \param mSq Invariant mass squared
  /// \param mR Mass of the resonant state
  /// \param ma Mass of daughter particle
  /// \param mb Mass of daughter particle
  /// \param width Decay width
  /// \param J Spin
  /// \param mesonRadius Meson Radius
  /// \param ffType Form factor type
  /// \return Amplitude value
  static std::complex<double>
  dynamicalFunction(double mSq, double mR, double ma, double mb, double width,
                    unsigned int J, double mesonRadius, formFactorType ffType);

  //============ SET/GET =================

  void SetWidthParameter(std::shared_ptr<ComPWA::FitParameter> w) {
    _width = w;
  }

  std::shared_ptr<ComPWA::FitParameter> GetWidthParameter() {
    return _width;
  }

  void SetWidth(double w) { _width->setValue(w); }

  double GetWidth() const { return _width->value(); }

  void SetMesonRadiusParameter(std::shared_ptr<ComPWA::FitParameter> r) {
    _mesonRadius = r;
  }

  std::shared_ptr<ComPWA::FitParameter> GetMesonRadiusParameter() {
    return _mesonRadius;
  }

  /// \see GetMesonRadius() const { return _mesonRadius->value(); }
  void SetMesonRadius(double w) { _mesonRadius->setValue(w); }

  /// Get meson radius.
  /// The meson radius is a measure of the size of the resonant state. It is
  /// used to calculate the angular momentum barrier factors.
  double GetMesonRadius() const { return _mesonRadius->value(); }

  /// \see GetFormFactorType()
  void SetFormFactorType(formFactorType t) { _ffType = t; }

  /// Get form factor type.
  /// The type of formfactor that is used to calculate the angular momentum
  /// barrier factors.
  formFactorType GetFormFactorType() { return _ffType; }

  virtual void GetParameters(ComPWA::ParameterList &list);

  virtual void GetParametersFast(std::vector<double> &list) const {
    AbstractDynamicalFunction::GetParametersFast(list);
    list.push_back(GetWidth());
    list.push_back(GetMesonRadius());
  }

  /// Update parameters to the values given in \p par
  virtual void updateParameters(const ComPWA::ParameterList &par);

  //=========== FUNCTIONTREE =================

  virtual bool hasTree() const { return true; }

  virtual std::shared_ptr<ComPWA::FunctionTree>
  tree(const ParameterList &sample, int pos, std::string suffix = "");

protected:
  /// Decay width of resonante state
  std::shared_ptr<ComPWA::FitParameter> _width;

  /// Meson radius of resonant state
  std::shared_ptr<ComPWA::FitParameter> _mesonRadius;

  /// Form factor type
  formFactorType _ffType;

private:
  /// Temporary values (used to trigger recalculation of normalization)
  double _current_mesonRadius;
  double _current_width;
};

class BreitWignerStrategy : public ComPWA::Strategy {
public:
  BreitWignerStrategy(std::string namee = "")
      : ComPWA::Strategy(ParType::MCOMPLEX), name(namee) {}

  virtual const std::string to_str() const {
    return ("relativistic BreitWigner of " + name);
  }

  virtual void execute(ComPWA::ParameterList &paras,
                       std::shared_ptr<ComPWA::Parameter> &out);

protected:
  std::string name;
};

} // namespace DecayDynamics
} // namespace Physics
} // namespace ComPWA

#endif
