// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_COEFFICIENTAMPLITUDEDECORATOR_HPP_
#define COMPWA_PHYSICS_COEFFICIENTAMPLITUDEDECORATOR_HPP_

#include "Amplitude.hpp"

namespace ComPWA {
namespace Physics {

class CoefficientAmplitudeDecorator : public NamedAmplitude {
public:
  CoefficientAmplitudeDecorator(
      const std::string &name, std::shared_ptr<Amplitude> Amplitude_,
      std::shared_ptr<ComPWA::FunctionTree::FitParameter> Magnitude_,
      std::shared_ptr<ComPWA::FunctionTree::FitParameter> Phase_)
      : NamedAmplitude(name), UndecoratedAmplitude(Amplitude_),
        Magnitude(Magnitude_), Phase(Phase_) {}

  std::complex<double> evaluate(const DataPoint &point) const final {
    return std::polar(Magnitude->value(), Phase->value()) *
           UndecoratedAmplitude->evaluate(point);
  }

  void addUniqueParametersTo(ComPWA::FunctionTree::ParameterList &list) final {
    Magnitude = list.addUniqueParameter(Magnitude);
    Phase = list.addUniqueParameter(Phase);
    UndecoratedAmplitude->addUniqueParametersTo(list);
  }

  void addFitParametersTo(std::vector<double> &FitParameters) final {
    FitParameters.push_back(Magnitude->value());
    FitParameters.push_back(Phase->value());
    UndecoratedAmplitude->addFitParametersTo(FitParameters);
  }

  void
  updateParametersFrom(const ComPWA::FunctionTree::ParameterList &list) final {
    std::shared_ptr<ComPWA::FunctionTree::FitParameter> mag =
        FindParameter(Magnitude->name(), list);
    Magnitude->updateParameter(mag);
    std::shared_ptr<ComPWA::FunctionTree::FitParameter> phase =
        FindParameter(Phase->name(), list);
    Phase->updateParameter(phase);

    UndecoratedAmplitude->updateParametersFrom(list);
  }

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createFunctionTree(const ComPWA::FunctionTree::ParameterList &DataSample,
                     const std::string &suffix) const final {

    size_t n = DataSample.mDoubleValue(0)->values().size();

    std::string nodeName = "CoefficientAmplitude(" + getName() + ")" + suffix;

    auto tr = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
        nodeName, ComPWA::FunctionTree::MComplex("", n),
        std::make_shared<ComPWA::FunctionTree::MultAll>(
            ComPWA::FunctionTree::ParType::MCOMPLEX));
    tr->createNode(
        "Strength",
        std::make_shared<ComPWA::FunctionTree::Value<std::complex<double>>>(),
        std::make_shared<ComPWA::FunctionTree::Complexify>(
            ComPWA::FunctionTree::ParType::COMPLEX),
        nodeName);
    tr->createLeaf("Magnitude", Magnitude, "Strength");
    tr->createLeaf("Phase", Phase, "Strength");

    tr->insertTree(UndecoratedAmplitude->createFunctionTree(DataSample, suffix),
                   nodeName);
    return tr;
  }

private:
  std::shared_ptr<Amplitude> UndecoratedAmplitude;
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> Magnitude;
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> Phase;
};

} // namespace Physics
} // namespace ComPWA

#endif
