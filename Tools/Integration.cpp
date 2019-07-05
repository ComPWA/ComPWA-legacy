// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>

#include "Core/FunctionTree/Intensity.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Logging.hpp"
#include "Data/DataSet.hpp"
#include "Integration.hpp"

#include "ThirdParty/parallelstl/include/pstl/algorithm"
#include "ThirdParty/parallelstl/include/pstl/execution"

namespace ComPWA {
namespace Tools {

/*
struct testGauss {
  double mu = 0;
  double sigma = 1;
  double operator()(double x) const {
    // Normalized Gaussian
    double expon = (-0.5) * (mu - x) * (mu - x) / (sigma * sigma);
    double val = 1 / (sigma * std::sqrt(2 * M_PI)) * std::exp(expon);
    return val;
  }
};

template <typename T> class IntegralByQuadrature {

public:
  IntegralByQuadrature(const T &func, std::pair<double, double> lim)
      : _func(func), _limits(lim), _depth(1), _integral(0.0) {}

  double Integral(int precision = 100) {
    while (_depth < precision) {
      //          std::cout<< "Current integral approximation: "<<_integral
      //          <<" depth: "<<_depth<<std::endl;
      Next();
      _depth *= 2;
    }
    return _integral;
  }

protected:
  const T &_func;
  std::pair<double, double> _limits;
  int _depth;
  double _integral;

  double Next() {
    double range = (_limits.second - _limits.first);
    if (_depth == 1) {
      _integral =
          0.5 * (range) * (_func(_limits.first) + _func(_limits.second));
    } else {
      double s = 0;
      int n = 0;
      double stepSize = 2 * (_limits.second - _limits.first) / (_depth);
      double x = _limits.first + 0.5 * stepSize;
      while (x < _limits.second) {
        s += _func(x);
        x += stepSize;
        n++;
      }
      if (n <= 0)
        throw std::runtime_error(
            "Tools::IntegralByQuadrature::Next() | Dividsion by zero!");
      _integral = 0.5 * (_integral + range * s / n);
    }
    return _integral;
  }
};*/

MCIntegrationStrategy::MCIntegrationStrategy(
    ComPWA::FunctionTree::ParameterList PhspDataSampleList_, double phspvolume)
    : PhspDataSampleList(PhspDataSampleList_), PhspVolume(phspvolume) {}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
MCIntegrationStrategy::createFunctionTree(
    std::shared_ptr<const ComPWA::FunctionTree::OldIntensity> intensity,
    const std::string &suffix) const {

  using namespace ComPWA::FunctionTree;
  double PhspWeightSum(PhspDataSampleList.mDoubleValue(0)->values().size());

  std::shared_ptr<Value<std::vector<double>>> phspweights;
  try {
    phspweights = findMDoubleValue("Weight", PhspDataSampleList);
    PhspWeightSum = std::accumulate(phspweights->values().begin(),
                                    phspweights->values().end(), 0.0);
  } catch (const Exception &e) {
  }

  auto tr = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      "Normalization", ValueFactory(ParType::DOUBLE),
      std::shared_ptr<Strategy>(new Inverse(ParType::DOUBLE)));
  tr->createNode("Integral",
                 std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE)),
                 "Normalization");
  // normTree->createLeaf("PhspVolume", PhspVolume, "Integral");
  tr->createLeaf("InverseSampleWeights", 1.0 / PhspWeightSum, "Integral");
  tr->createNode("Sum", std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)),
                 "Integral");
  tr->createNode("WeightedIntensities",
                 std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)),
                 "Sum");

  if (phspweights)
    tr->createLeaf("EventWeight", phspweights, "WeightedIntensities");
  tr->insertTree(intensity->createFunctionTree(PhspDataSampleList, suffix),
                 "WeightedIntensities");

  return tr;
}

double integrate(std::shared_ptr<Intensity> intensity,
                 ComPWA::Data::DataSet phspsample, double phspVolume) {

  std::vector<double> Intensities = intensity->evaluate(phspsample.Data);
  std::transform(
      pstl::execution::par_unseq, Intensities.begin(), Intensities.end(),
      phspsample.Weights.begin(), Intensities.begin(),
      [](double intensity, double weight) { return intensity * weight; });
  double IntensitySum(
      std::accumulate(Intensities.begin(), Intensities.end(), 0.0));
  double WeightSum(std::accumulate(phspsample.Weights.begin(),
                                   phspsample.Weights.end(), 0.0));

  return (IntensitySum * phspVolume / WeightSum);
}

double maximum(std::shared_ptr<Intensity> intensity,
               ComPWA::Data::DataSet Sample) {
  if (!Sample.Weights.size()) {
    LOG(DEBUG) << "Tools::Maximum(): Maximum can not be determined since "
                  "sample is empty.";
    return 1.0;
  }

  std::vector<double> Intensities = intensity->evaluate(Sample.Data);

  std::transform(Intensities.begin(), Intensities.end(), Sample.Weights.begin(),
                 Intensities.begin(), [](double Intensity, double Weight) {
                   return Intensity * Weight;
                 });
  // determine maximum
  double max(*std::max_element(Intensities.begin(), Intensities.end()));
  LOG(DEBUG) << "Tools::Maximum(): found maximum value of " << max;
  return max;
}

} // namespace Tools
} // namespace ComPWA
