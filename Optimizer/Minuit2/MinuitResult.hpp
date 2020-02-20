// Copyright (c) 2014, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef OPTIMIZER_MINUIT2_MINUITRESULT_HPP_
#define OPTIMIZER_MINUIT2_MINUITRESULT_HPP_

#include "Core/FitResult.hpp"

#include "boost/serialization/base_object.hpp"

namespace ROOT {
namespace Minuit2 {
class FunctionMinimum;
} // namespace Minuit2
} // namespace ROOT

namespace ComPWA {
namespace Optimizer {
namespace Minuit2 {

struct MinuitResult : public FitResult {
  MinuitResult() = default;
  MinuitResult(const FitResult &Result,
               const ROOT::Minuit2::FunctionMinimum &FMin);
  bool CovPosDef = false;           // covariance matrix pos.-def.
  bool HasValidParameters = false;  // valid parameters
  bool HasValidCov = false;         // valid covariance
  bool HasAccCov = false;           // accurate covariance
  bool HasReachedCallLimit = false; // call limit reached
  bool EdmAboveMax = false;
  bool HesseFailed = false;
  double ErrorDef = false;
  unsigned int NFcn = 0;
  double Edm = 0.0; // estimated distance to minimum
  std::vector<double> GlobalCC;

  void write(std::string filename) const;

  friend std::ostream &operator<<(std::ostream &os, const MinuitResult &Result);

private:
  friend class boost::serialization::access;
  template <class archive>
  void serialize(archive &ar, const unsigned int version) {
    using namespace boost::serialization;
    ar &BOOST_SERIALIZATION_BASE_OBJECT_NVP(FitResult);
    ar &BOOST_SERIALIZATION_NVP(CovPosDef);
    ar &BOOST_SERIALIZATION_NVP(HasValidParameters);
    ar &BOOST_SERIALIZATION_NVP(HasValidCov);
    ar &BOOST_SERIALIZATION_NVP(HasAccCov);
    ar &BOOST_SERIALIZATION_NVP(HasReachedCallLimit);
    ar &BOOST_SERIALIZATION_NVP(EdmAboveMax);
    ar &BOOST_SERIALIZATION_NVP(HesseFailed);
    ar &BOOST_SERIALIZATION_NVP(ErrorDef);
    ar &BOOST_SERIALIZATION_NVP(NFcn);
    ar &BOOST_SERIALIZATION_NVP(Edm);
    ar &BOOST_SERIALIZATION_NVP(GlobalCC);
  }
};

MinuitResult load(std::string filename);

} // namespace Minuit2
} // namespace Optimizer
} // namespace ComPWA

#endif
