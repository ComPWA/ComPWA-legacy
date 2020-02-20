// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef CORE_FITRESULT_HPP_
#define CORE_FITRESULT_HPP_

#include "Core/FitParameter.hpp"
#include "Core/Function.hpp"

#include <boost/serialization/vector.hpp>

#include <chrono>

namespace ComPWA {

/// Data structure which resembles a general fit result. Optimizers should
/// derive from this structure and append more information via inheritance.
struct FitResult {
  FitParameterList InitialParameters;
  FitParameterList FinalParameters;
  unsigned int NumFreeParameters;
  bool IsValid = false;

  double InitialEstimatorValue = 0.0;
  double FinalEstimatorValue = 0.0;
  std::chrono::seconds FitDuration = std::chrono::seconds(0);

  std::vector<std::vector<double>> CovarianceMatrix;

  friend std::ostream &operator<<(std::ostream &os, const FitResult &Result);

  void write(std::string filename) const;

private:
  friend class boost::serialization::access;
  template <class archive>
  void serialize(archive &ar, const unsigned int version) {
    ar &BOOST_SERIALIZATION_NVP(InitialParameters);
    ar &BOOST_SERIALIZATION_NVP(FinalParameters);
    ar &BOOST_SERIALIZATION_NVP(NumFreeParameters);
    ar &BOOST_SERIALIZATION_NVP(IsValid);
    ar &BOOST_SERIALIZATION_NVP(InitialEstimatorValue);
    ar &BOOST_SERIALIZATION_NVP(FinalEstimatorValue);
    auto x = FitDuration.count();
    ar &boost::serialization::make_nvp("FitDurationInSeconds", x);
    FitDuration = std::chrono::seconds(x);
    ar &BOOST_SERIALIZATION_NVP(CovarianceMatrix);
  }
};

FitResult load(std::string filename);

std::string makeFitParameterString(ComPWA::FitParameter<double> p);

void initializeWithFitResult(ComPWA::Intensity &Intens,
                             ComPWA::FitResult Result);
} // namespace ComPWA

#endif
