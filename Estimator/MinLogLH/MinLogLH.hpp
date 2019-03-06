// Copyright (c) 2013, 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_ESTIMATOR_MINLOGLH_HPP_
#define COMPWA_ESTIMATOR_MINLOGLH_HPP_

#include <memory>
#include <vector>

#include "Estimator/Estimator.hpp"
#include "Estimator/FunctionTreeEstimator.hpp"

namespace ComPWA {

class Intensity;
struct DataPoint;

namespace Data {
class DataSet;
}

namespace Estimator {

///
/// \class MinLogLH
/// Negative Log Likelihood-Estimator. This class calculates the negative log
/// likelihood -log(LH) using an Intensity and data samples. Data data samples
/// are retrieved from the DataStorage.
///
/// \par log likelihood
/// The negative log LH is given by:
/// \f[
///    -log \mathcal{L} = N_{\mathrm{obs}} * \log(\lambda) -
///    \sum_i^{N_{\mathrm{obs}}} \log(I(x_i))
/// \f]
/// with the Intensity \f$I(x_i)\f$ for a given event \f$x_i\f$ and the
/// phase integral
/// \f$ \lambda = \frac{V}{\sum w_i}\sum_i^{N_{\mathrm{phsp}}} I(x_i)\cdot
/// \epsilon_i \cdot w_i \f$ \f$ \epsilon \f$ and \f$ w \f$ are the efficiency
/// and weight for each event \f$ V \f$ is the phase space volume (in which the
/// phase space events lie)
///
/// The sum over all weights is necessary to normalize the weights to one.
/// Otherwise the error estimate would be incorrect. The Intensity does not
/// have to be normalized. This is done automatically by the phase space
/// integral.
///
/// \par Efficiency correction
/// It is assumed that the data already includes the efficiency.
///
class MinLogLH : public ComPWA::Estimator::Estimator {

public:
  MinLogLH(std::shared_ptr<ComPWA::Intensity> intensity,
           const std::vector<ComPWA::DataPoint> &datapoints,
           const std::vector<ComPWA::DataPoint> &phsppoints);

  /// Value of log likelihood function.
  double evaluate() const final;

private:
  std::shared_ptr<ComPWA::Intensity> Intensity;

  const std::vector<DataPoint> &DataPoints;
  const std::vector<DataPoint> &PhspDataPoints;
};

std::shared_ptr<FunctionTree> createMinLogLHEstimatorFunctionTree(
    std::shared_ptr<ComPWA::Intensity> Intensity,
    std::shared_ptr<ComPWA::Data::DataSet> DataSample,
    std::shared_ptr<ComPWA::Data::DataSet> PhspDataSample = {});

std::shared_ptr<FunctionTreeEstimator> createMinLogLHFunctionTreeEstimator(
    std::shared_ptr<ComPWA::Intensity> Intensity,
    std::shared_ptr<ComPWA::Data::DataSet> DataSample,
    std::shared_ptr<ComPWA::Data::DataSet> PhspDataSample = {});

} // namespace Estimator
} // namespace ComPWA

#endif
