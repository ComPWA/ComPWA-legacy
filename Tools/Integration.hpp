// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_TOOLS_INTEGRATION_HPP_
#define COMPWA_TOOLS_INTEGRATION_HPP_

#include <memory>
#include <vector>

#include "Core/FunctionTree.hpp"
#include "Core/Function.hpp"

namespace ComPWA {

class OldIntensity;
struct DataPoint;
class Kinematics;

namespace Data {
class DataSet;
}

namespace Tools {

class IntegrationStrategy {
public:
  virtual ~IntegrationStrategy() = default;

  virtual double
  integrate(std::shared_ptr<ComPWA::OldIntensity> intensity) const = 0;

  virtual std::shared_ptr<ComPWA::FunctionTree>
  createFunctionTree(std::shared_ptr<const ComPWA::OldIntensity> intensity,
                     const std::string &suffix) const = 0;
};

class MCIntegrationStrategy : public IntegrationStrategy {
public:
  MCIntegrationStrategy(std::shared_ptr<const ComPWA::Data::DataSet> phspsample,
                        double phspvolume = 1.0);
  double
  integrate(std::shared_ptr<ComPWA::OldIntensity> intensity) const final;

  std::shared_ptr<ComPWA::FunctionTree>
  createFunctionTree(std::shared_ptr<const ComPWA::OldIntensity> intensity,
                     const std::string &suffix) const final;

private:
  std::shared_ptr<const ComPWA::Data::DataSet> PhspSample;
  double PhspVolume;
};

double integrate(std::shared_ptr<OldIntensity> intensity,
                 std::shared_ptr<const ComPWA::Data::DataSet> phspsample,
                 double phspVolume = 1.0);

double maximum(std::shared_ptr<Intensity> intensity,
               const std::vector<DataPoint> &sample);

double maximum(std::shared_ptr<Intensity> intensity,
               std::shared_ptr<const Data::DataSet> sample,
               std::shared_ptr<Kinematics> kin);

} // namespace Tools
} // namespace ComPWA

#endif
