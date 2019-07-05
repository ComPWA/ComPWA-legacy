// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_TOOLS_INTEGRATION_HPP_
#define COMPWA_TOOLS_INTEGRATION_HPP_

#include <memory>
#include <vector>

#include "Core/Function.hpp"
#include "Core/FunctionTree/FunctionTree.hpp"

namespace ComPWA {

class OldIntensity;
struct DataPoint;
class Kinematics;

namespace Data {
struct DataSet;
}

namespace Tools {

class IntegrationStrategy {
public:
  virtual ~IntegrationStrategy() = default;

  virtual std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createFunctionTree(
      std::shared_ptr<const ComPWA::FunctionTree::OldIntensity> intensity,
      const std::string &suffix) const = 0;
};

class MCIntegrationStrategy : public IntegrationStrategy {
public:
  MCIntegrationStrategy(ComPWA::FunctionTree::ParameterList PhspDataSampleList_,
                        double phspvolume = 1.0);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree> createFunctionTree(
      std::shared_ptr<const ComPWA::FunctionTree::OldIntensity> intensity,
      const std::string &suffix) const final;

private:
  ComPWA::FunctionTree::ParameterList PhspDataSampleList;
  double PhspVolume;
};

double integrate(std::shared_ptr<Intensity> intensity,
                 ComPWA::Data::DataSet phspsample, double phspVolume = 1.0);

double maximum(std::shared_ptr<Intensity> intensity,
               ComPWA::FunctionTree::ParameterList sample);

double maximum(std::shared_ptr<Intensity> intensity,
               std::shared_ptr<const Data::DataSet> sample,
               std::shared_ptr<Kinematics> kin);

} // namespace Tools
} // namespace ComPWA

#endif
