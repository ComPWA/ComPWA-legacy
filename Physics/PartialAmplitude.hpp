// Copyright (c) 2016, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Contains PartialAmplitude interface class
///

#ifndef PHYSICS_PartialAmplitude_HPP_
#define PHYSICS_PartialAmplitude_HPP_

#include <vector>
#include <memory>

#include "Core/FitParameter.hpp"
#include "Core/FitParameter.hpp"
#include "Core/DataPoint.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Kinematics.hpp"

namespace ComPWA {
namespace Physics {

///
/// \class PartialAmplitude
/// PartialAmplitude class is an interface classe that resembles a two-body
/// decay.
///
class PartialAmplitude {

public:
  PartialAmplitude()
      : PreFactor(1, 0), PhspVolume(1), CurrentIntegral(1.0),
        CurrentMagnitude(0.0), CurrentPhase(0.0){};

  virtual ~PartialAmplitude(){};

  virtual PartialAmplitude *clone(std::string newName = "") const = 0;

  virtual boost::property_tree::ptree save() const = 0;

  /// Get current normalization.
  virtual double normalization() const = 0;

  /// Check of parameters have changed and normalization has to be recalculatecd
  virtual bool isModified() const = 0;

  /// Label as modified/unmodified
  virtual void setModified(bool b) = 0;

  /// Value of PartialAmplitude at \param point with normalization factor
  virtual std::complex<double> evaluate(const DataPoint &point) const {
    return evaluateNoNorm(point) * normalization() * coefficient() * prefactor();
  }

  /// Value of PartialAmplitude at \param point without normalization factor
  virtual std::complex<double> evaluateNoNorm(const DataPoint &point) const = 0;

  virtual std::string name() const { return Name; }

  virtual void setName(std::string name) { Name = name; }

  virtual void setPrefactor(std::complex<double> pre) { PreFactor = pre; }

  virtual std::complex<double> prefactor() const { return PreFactor; }

  virtual std::complex<double> coefficient() const {
    return std::polar(magnitude(), phase());
  }

  std::shared_ptr<ComPWA::FitParameter> magnitudeParameter() {
    return Magnitude;
  }

  double magnitude() const { return std::fabs(Magnitude->value()); }

  void setMagnitudeParameter(std::shared_ptr<ComPWA::FitParameter> par) {
    Magnitude = par;
  }

  void setMagnitude(double par) { Magnitude->setValue(par); }

  std::shared_ptr<ComPWA::FitParameter> phaseParameter() { return Phase; }

  double phase() const { return Phase->value(); }

  void setPhaseParameter(std::shared_ptr<ComPWA::FitParameter> par) {
    Phase = par;
  }

  void setPhase(double par) { Phase->setValue(par); }

  virtual void parameters(ParameterList &list) {
    Phase = list.addUniqueParameter(Phase);
    Magnitude = list.addUniqueParameter(Magnitude);
  }

  /// Fill vector with parameters (fast). No check is performed if parameters
  /// already exist. \see parameters(ParameterList &list)
  virtual void parametersFast(std::vector<double> &list) const {
    list.push_back(magnitude());
    list.push_back(phase());
  }

  /// Update parameters to the values given in \p list
  virtual void updateParameters(const ParameterList &list) {
    // Try to update magnitude
    std::shared_ptr<FitParameter> mag;
    try {
      mag = FindParameter(Magnitude->name(), list);
    } catch (std::exception &ex) {
    }
    if (mag)
      Magnitude->updateParameter(mag);
    std::shared_ptr<FitParameter> phase;

    // Try to update phase
    try {
      phase = FindParameter(Phase->name(), list);
    } catch (std::exception &ex) {
    }
    if (phase)
      Phase->updateParameter(phase);
  }

  /// Set phase space sample.
  /// We use the phase space sample to calculate the normalization. The sample
  /// should be without efficiency applied.
  virtual void
  setPhspSample(std::shared_ptr<std::vector<ComPWA::DataPoint>> phspSample) {
    PhspSample = phspSample;
  };

  virtual void setPhspVolume(double phspVol) { PhspVolume = phspVol; }

  virtual double phspVolume() const { return PhspVolume; }

  virtual bool hasTree() const { return false; }

  virtual std::shared_ptr<FunctionTree>
  tree(std::shared_ptr<Kinematics> kin, const ComPWA::ParameterList &sample,
       const ComPWA::ParameterList &toySample, std::string suffix) = 0;

protected:
  std::string Name;
  std::shared_ptr<ComPWA::FitParameter> Magnitude;
  std::shared_ptr<ComPWA::FitParameter> Phase;
  std::complex<double> PreFactor;

  /// Phsp sample for numerical integration
  std::shared_ptr<std::vector<ComPWA::DataPoint>> PhspSample;

  /// The phase-space volume is needed for proper normalization of the
  /// PartialAmplitude
  double PhspVolume;

  virtual double integral() const {
    if (!PhspSample->size()) {
      LOG(DEBUG)
          << "PartialAmplitude::Integral() | Integral can not be "
             " calculated since no phsp sample is set. "
             " Set a sample using SetPhspSamples(phspSample, toySample)!";
      return 1.0;
    }

    double sumIntens = 0;
    for (auto i : *PhspSample.get()) {
      sumIntens += std::norm(evaluateNoNorm(i));
    }

    double integral = (sumIntens * PhspVolume / PhspSample->size());
    LOG(TRACE) << "PartialAmplitude::Integral() | Integral is " << integral
               << ".";
    assert(!std::isnan(integral));
    return integral;
  }

  /// Integral value (temporary)
  double CurrentIntegral;

  /// Temporary value of magnitude (used to trigger recalculation of
  /// normalization)
  double CurrentMagnitude;
  /// Temporary value of phase (used to trigger recalculation of normalization)
  double CurrentPhase;
};

///
/// \class NonResonant
/// Decay to a multi particle final state without resonance dynamics.
///
class NonResonant : public ComPWA::Physics::PartialAmplitude {
public:
  NonResonant(std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
              const boost::property_tree::ptree &pt) {
    load(partL, kin, pt);
  }

  virtual NonResonant *clone(std::string newName = "") const {
    auto tmp = new NonResonant(*this);
    tmp->setName(newName);
    return tmp;
  }

  void load(std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
            const boost::property_tree::ptree &pt) {

    LOG(TRACE) << "NonResonant::load() |";
    setPhspVolume(kin->phspVolume());

    Name = pt.get<std::string>("<xmlattr>.Name", "empty");
    Magnitude =
        std::make_shared<ComPWA::FitParameter>("Magnitude_" + Name, 1.0);
    Phase = std::make_shared<ComPWA::FitParameter>("Phase_" + Name, 0.0);
    std::shared_ptr<FitParameter> mag, phase;
    for (const auto &v : pt.get_child("")) {
      if (v.first == "Parameter") {
        if (v.second.get<std::string>("<xmlattr>.Type") == "Magnitude") {
          Magnitude = std::make_shared<FitParameter>(v.second);
        }
        if (v.second.get<std::string>("<xmlattr>.Type") == "Phase") {
          Phase = std::make_shared<FitParameter>(v.second);
        }
      } else {
        // ignored further settings. Should we throw an error?
      }
    }
  }

  virtual boost::property_tree::ptree save() const {
    boost::property_tree::ptree pt;
    pt.put<std::string>("<xmlattr>.Name", name());

    boost::property_tree::ptree tmp = Magnitude->save();
    tmp.put("<xmlattr>.Type", "Magnitude");
    pt.add_child("Parameter", tmp);

    tmp = Phase->save();
    tmp.put("<xmlattr>.Type", "Phase");
    pt.add_child("Parameter", tmp);

    return pt;
  }

  /// Check of parameters have changed and normalization has to be recalculatecd
  virtual bool isModified() const {
    if (magnitude() != CurrentMagnitude || phase() != CurrentPhase)
      return true;
    return false;
  }

  /// Label as modified/unmodified
  virtual void setModified(bool b) {
    if (b) {
      CurrentMagnitude = std::numeric_limits<double>::quiet_NaN();
      CurrentPhase = std::numeric_limits<double>::quiet_NaN();
    } else {
      CurrentMagnitude = magnitude();
      CurrentPhase = phase();
    }
  }

  /// Value of PartialAmplitude at \param point without normalization factor
  virtual std::complex<double> evaluateNoNorm(const DataPoint &point) const {
    return std::complex<double>(1., 0.);
  };

  /// Get current normalization.
  virtual double normalization() const { return 1 / std::sqrt(PhspVolume); };

  //=========== FUNCTIONTREE =================
  virtual bool hasTree() const { return true; }

  virtual std::shared_ptr<FunctionTree>
  tree(std::shared_ptr<Kinematics> kin, const ComPWA::ParameterList &sample,
       const ComPWA::ParameterList &toySample, std::string suffix) {

    size_t n = sample.mDoubleValue(0)->values().size();

    std::string nodeName = "PartialAmplitude(" + name() + ")" + suffix;

    auto tr = std::make_shared<FunctionTree>(
        nodeName, MComplex("", n),
        std::make_shared<MultAll>(ParType::MCOMPLEX));
    tr->createNode("Strength", std::make_shared<Value<std::complex<double>>>(),
                   std::make_shared<Complexify>(ParType::COMPLEX), nodeName);
    tr->createLeaf("Magnitude", Magnitude, "Strength");
    tr->createLeaf("Phase", Phase, "Strength");
    tr->createLeaf("PreFactor", PreFactor, nodeName);
    tr->createLeaf("One", MComplex("", n, std::complex<double>(1., 0.)),
                   nodeName);

    tr->createNode("Normalization",
                   std::make_shared<Value<double>>(1 / std::sqrt(PhspVolume)),
                   std::make_shared<Inverse>(ParType::DOUBLE),
                   nodeName); // 1/normLH

    tr->parameter();
    return tr;
  };
};

} // ns::Physics
} // ns::ComPWA
#endif
