// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Minuit2 interface FCN base class.
///

#ifndef COMPWA_OPTIMIZER_MINUIT2_MINUITFCN_HPP_
#define COMPWA_OPTIMIZER_MINUIT2_MINUITFCN_HPP_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Core/ParameterList.hpp"

#include "Minuit2/FCNBase.h"

namespace ComPWA {
namespace Estimator {
class Estimator;
}
} // namespace ComPWA

namespace ROOT {
namespace Minuit2 {

///
/// \class MinuitFcn
/// Minuit2 function to be optimized based on the Minuit2 FcnBase. This class
/// uses the Estimator interface for the optimization.
///
class MinuitFcn : public FCNBase {

public:
  MinuitFcn(std::shared_ptr<ComPWA::Estimator::Estimator> estimator,
            ComPWA::ParameterList &parameters);
  virtual ~MinuitFcn();

  double operator()(const std::vector<double> &x) const;

  double Up() const;

  inline void setNameID(const unsigned int id, const std::string &name) {
    auto result = IDToParameterNameMapping.insert(
        std::pair<unsigned int, std::string>(id, name));
    if (!result.second) {
      std::stringstream ss;
      ss << "MinuitFcn::setNameID(): Could not create entry in ID-name map for "
            "id="
         << id << " and name=" << name;
      throw std::runtime_error(ss.str());
    }
  };

  inline std::string parName(const unsigned int id) {
    return IDToParameterNameMapping.at(id);
  };

private:
  std::shared_ptr<ComPWA::Estimator::Estimator> Estimator;

  /// List of Parameters that influence the Estimator
  ComPWA::ParameterList &Parameters;

  /// mapping of minuit ids to ComPWA names
  std::map<unsigned int, std::string> IDToParameterNameMapping;
};

} // namespace Minuit2
} // namespace ROOT

#endif
