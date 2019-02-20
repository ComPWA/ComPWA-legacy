// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_TOOLS_INTEGRATION_HPP_
#define COMPWA_TOOLS_INTEGRATION_HPP_

#include <memory>
#include <vector>

#include "Core/FunctionTree.hpp"

namespace ComPWA {

class Intensity;
class DataPoint;
class Kinematics;
class Event;

namespace Tools {

class IntegrationStrategy {
public:
  virtual ~IntegrationStrategy() = default;

  virtual double
  integrate(std::shared_ptr<const ComPWA::Intensity> intensity) const = 0;

  virtual std::shared_ptr<ComPWA::FunctionTree>
  createFunctionTree(std::shared_ptr<const ComPWA::Intensity> intensity,
                     const std::string &suffix) const = 0;
};

class MCIntegrationStrategy : public IntegrationStrategy {
public:
  MCIntegrationStrategy(const std::vector<DataPoint> &phsppoints,
                        double phspvolume = 1.0);
  double
  integrate(std::shared_ptr<const ComPWA::Intensity> intensity) const final;

  std::shared_ptr<ComPWA::FunctionTree>
  createFunctionTree(std::shared_ptr<const ComPWA::Intensity> intensity,
                     const std::string &suffix) const final;

private:
  const std::vector<ComPWA::DataPoint> PhspDataPoints;
  double PhspVolume;
};

double integrate(std::shared_ptr<const Intensity> intensity,
                 const std::vector<DataPoint> &sample, double phspVolume = 1.0);

double maximum(std::shared_ptr<const Intensity> intensity,
               const std::vector<DataPoint> &sample);

double maximum(std::shared_ptr<const Intensity> intensity,
               const std::vector<Event> &sample,
               std::shared_ptr<Kinematics> kin);

} // namespace Tools
} // namespace ComPWA

#endif
