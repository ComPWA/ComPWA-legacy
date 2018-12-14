// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_TOOLS_INTEGRATION_HPP_
#define COMPWA_TOOLS_INTEGRATION_HPP_

#include <cmath>
#include <complex>

#include "Core/Intensity.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Logging.hpp"
#include "Data/DataTransformation.hpp"

namespace ComPWA {
namespace Tools {

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
};

inline double Integral(std::shared_ptr<const Intensity> intens,
                       const std::vector<DataPoint> &sample,
                       double phspVolume = 1.0) {

  if (!sample.size()) {
    LOG(DEBUG) << "Tools::Integral() | Integral can not be calculated "
                  "since phsp sample is empty.";
    return 1.0;
  }
  double IntensitySum(0.0);
  for (auto i : sample)
    IntensitySum += intens->evaluate(i);

  return (IntensitySum * phspVolume / sample.size());
}

inline double Integral(std::shared_ptr<const Intensity> intens,
                       std::shared_ptr<std::vector<DataPoint>> sample,
                       double phspVolume = 1.0) {
  return Integral(intens, *sample.get(), phspVolume);
}

inline double Maximum(std::shared_ptr<Intensity> intens,
                      std::shared_ptr<std::vector<DataPoint>> sample) {

  if (!sample->size()) {
    LOG(DEBUG) << "Tools::Maximum() | Maximum can not be determined since "
                  "sample is empty.";
    return 1.0;
  }

  double max = 0;
  for (auto i : *sample.get()) {
    double val = intens->evaluate(i);
    if (val > max)
      max = val;
  }

  return max;
}

inline double Maximum(std::shared_ptr<Kinematics> kin,
                      std::shared_ptr<Intensity> intens,
                      const std::vector<Event> &sample) {

  if (0 == sample.size()) {
    LOG(DEBUG)
        << "Maximum() | MAximum can not be determined since sample is empty.";
    return 0.0;
  }

  auto DataPoints = Data::convertEventsToDataPoints(sample, kin);

  double max(0.0);
  DataPoint maxPoint;
  for (auto const &x : DataPoints) {
    double val = intens->evaluate(x);
    if (val > max) {
      maxPoint = x;
      max = val;
    }
  }

  LOG(DEBUG) << "Maximum() | Maximum found at " << maxPoint << ".";
  return max;
}

} // namespace Tools
} // namespace ComPWA
#endif
