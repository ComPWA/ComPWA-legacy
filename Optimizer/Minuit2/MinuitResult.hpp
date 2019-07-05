// Copyright (c) 2014, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef OPTIMIZER_MINUIT2_MINUITRESULT_HPP_
#define OPTIMIZER_MINUIT2_MINUITRESULT_HPP_

#include "Core/FitResult.hpp"

#include "boost/serialization/base_object.hpp"

namespace ComPWA {
namespace Optimizer {
namespace Minuit2 {

struct MinuitResult : public FitResult {
  bool CalcInterference;
  bool IsValid;             // result valid
  bool CovPosDef;           // covariance matrix pos.-def.
  bool HasValidParameters;  // valid parameters
  bool HasValidCov;         // valid covariance
  bool HasAccCov;           // accurate covariance
  bool HasReachedCallLimit; // call limit reached
  bool EdmAboveMax;
  bool HesseFailed;
  double ErrorDef;
  unsigned int NFcn;
  double Edm; // estimated distance to minimum
  std::vector<double> GlobalCC;

  void print(std::ostream &os) const;
  friend std::ostream &operator<<(std::ostream &os, const MinuitResult &Result);

private:
  friend class boost::serialization::access;
  template <class archive>
  void serialize(archive &ar, const unsigned int version) {
    using namespace boost::serialization;
    ar &BOOST_SERIALIZATION_BASE_OBJECT_NVP(FitResult);
    ar &BOOST_SERIALIZATION_NVP(CalcInterference);
    ar &BOOST_SERIALIZATION_NVP(IsValid);
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

} // namespace Minuit2
} // namespace Optimizer
} // namespace ComPWA

#endif
