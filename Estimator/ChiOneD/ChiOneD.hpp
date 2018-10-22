// Copyright (c) 2013 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//! Simple \f$\chi^{2}\f$-Estimator.
/*! \class ChiOneD
 * @file ChiOneD.hpp
 * This class calculates a simple \f$\chi^{2}\f$ of a intensity and a dataset.
 * Data and Model are provided in the constructor using the Amplitude and Data
 * interfaces. The class itself fulfills the Estimator interface.
*/

#ifndef COMPWA_ESTIMATOR_CHIONED_HPP_
#define COMPWA_ESTIMATOR_CHIONED_HPP_
#include <vector>
#include <memory>
#include <string>

// PWA-Headers
#include "Core/Estimator.hpp"

namespace ComPWA {

namespace Data {
class Data;
}

class AmpIntensity;
class Event;
class ParameterList;
class FunctionTree;
class Kinematics;

namespace Estimator {
namespace ChiOneD {

class ChiOneD : public ComPWA::IEstimator {

public:
  ChiOneD(std::shared_ptr<Kinematics> kin, std::shared_ptr<AmpIntensity>,
          std::shared_ptr<Data::Data>);

  virtual double controlParameter(ParameterList &minPar);

  virtual bool hasTree() { return false; }
  
  virtual std::shared_ptr<FunctionTree> tree() {
    return std::shared_ptr<FunctionTree>();
  }
  
  virtual std::shared_ptr<AmpIntensity> GetIntensity() {
    return std::shared_ptr<AmpIntensity>();
  }

protected:

private:
  std::shared_ptr<Kinematics> kin_;
  std::shared_ptr<AmpIntensity> pPIF_;
  std::shared_ptr<Data::Data> pDIF_;
};

} /* namespace ChiOneD */
} /* namespace Estimator */
} /* namespace ComPWA */

#endif /* _EIFChiOneD_HPP */
