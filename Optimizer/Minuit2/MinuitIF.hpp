// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Minuit interface
///

#ifndef _MINUITIF_HPP
#define _MINUITIF_HPP

#include <vector>
#include <memory>

#include <boost/serialization/nvp.hpp>

#include "Minuit2/MnStrategy.h"

#include "Core/ParameterList.hpp"
#include "Optimizer/Optimizer.hpp"
#include "Optimizer/Minuit2/MinuitFcn.hpp"
#include "Optimizer/Minuit2/MinuitResult.hpp"

namespace ComPWA {
namespace Optimizer {
namespace Minuit2 {
///
/// \class MinuitIF
/// Wrapper of the Minuit2 Optimizer library. This class provides a wrapper
/// around the Minuit2 library. It fulfills the
/// Optimizer interface to be easily adapted to other modules. The data needs to
/// be provided with the ControlParameter interface.
///
class MinuitIF : public Optimizer {

public:
  MinuitIF(std::shared_ptr<ComPWA::IEstimator> esti, ParameterList& par);
  
  virtual std::shared_ptr<FitResult> exec(ParameterList& par);

  virtual ~MinuitIF();

  virtual void SetHesse(bool onoff) { enableHesse = onoff; }
  virtual bool GetHesse() { return enableHesse; }
  virtual void SetMinos(bool onoff) { enableMinos = onoff; }
  virtual bool GetMinos() { return enableMinos; }

protected:
private:
  ROOT::Minuit2::MinuitFcn _myFcn;
  std::shared_ptr<ComPWA::IEstimator> estimator;
  bool enableHesse;
  bool enableMinos;
};

class MinuitStrategy : public ROOT::Minuit2::MnStrategy {
public:
  MinuitStrategy(unsigned int i = 1) : MnStrategy(i) {
    fGradNCyc = GradientNCycles();
    fGradTlrStp = GradientStepTolerance();
    fGradTlr = GradientTolerance();
    fHessNCyc = HessianNCycles();
    fHessTlrStp = HessianStepTolerance();
    fHessTlrG2 = HessianG2Tolerance();
    fHessGradNCyc = HessianGradientNCycles();
  };
  void init() {
    SetGradientNCycles(fGradNCyc);
    SetGradientStepTolerance(fGradTlrStp);
    SetGradientTolerance(fGradTlr);
    SetHessianNCycles(fHessNCyc);
    SetHessianStepTolerance(fHessTlrStp);
    SetHessianG2Tolerance(fHessTlrG2);
    SetHessianGradientNCycles(fHessGradNCyc);
  }

private:
  unsigned int fGradNCyc;
  double fGradTlrStp;
  double fGradTlr;
  unsigned int fHessNCyc;
  double fHessTlrStp;
  double fHessTlrG2;
  unsigned int fHessGradNCyc;

  friend class boost::serialization::access;
  template <class archive>
  void serialize(archive &ar, const unsigned int version) {
    //    ar & BOOST_SERIALIZATION_NVP(fStrategy);
    ar &BOOST_SERIALIZATION_NVP(fGradNCyc);
    ar &BOOST_SERIALIZATION_NVP(fGradTlrStp);
    ar &BOOST_SERIALIZATION_NVP(fGradTlr);
    ar &BOOST_SERIALIZATION_NVP(fHessNCyc);
    ar &BOOST_SERIALIZATION_NVP(fHessTlrStp);
    ar &BOOST_SERIALIZATION_NVP(fHessTlrG2);
    ar &BOOST_SERIALIZATION_NVP(fHessGradNCyc);
  }
};

} // ns::Minuit2
} // ns::Optimizer
} // ns::ComPWA

#endif
