// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef OPTIMIZER_MINUIT2_MINUITIF_HPP
#define OPTIMIZER_MINUIT2_MINUITIF_HPP

#include "Core/FitParameter.hpp"
#include "Optimizer/Minuit2/MinuitFcn.hpp"
#include "Optimizer/Minuit2/MinuitResult.hpp"
#include "Optimizer/Optimizer.hpp"

namespace ROOT {
namespace Minuit2 {
class FunctionMinimum;
class MnStrategy;
} // namespace Minuit2
} // namespace ROOT

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
class MinuitIF : public Optimizer<MinuitResult> {
public:
  MinuitIF(bool UseHesse_ = true, bool UseMinos_ = false);

  void useMinos(bool x);
  void useHesse(bool x);

  MinuitResult optimize(ComPWA::Estimator::Estimator<double> &Estimator,
                        ComPWA::FitParameterList InitialParameters);

private:
  static FitParameterList
  getFinalParameters(const ROOT::Minuit2::MnUserParameterState &minState,
                     FitParameterList InitialParameters);
  static std::vector<std::vector<double>>
  getCovarianceMatrix(const ROOT::Minuit2::MnUserParameterState &minState);

  bool UseHesse;
  bool UseMinos;
};

} // namespace Minuit2
} // namespace Optimizer
} // namespace ComPWA

namespace boost {
namespace serialization {

template <class Archive>
void serialize(Archive &ar, ROOT::Minuit2::MnStrategy &s,
               const unsigned int version);
}
} // namespace boost

#endif
