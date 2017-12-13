// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Contains AmpWignerD class which provides WignerD functions.
///

#ifndef AMPWIGNER_D
#define AMPWIGNER_D

#include <vector>
#include <memory>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/exceptions.hpp>

#include "Core/ParameterList.hpp"
#include "Core/Functions.hpp"
#include "Core/DataPoint.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Spin.hpp"
#include "Core/Properties.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

///
/// \class AmpWignerD
/// Angular distribution based on WignerD functions
///
class AmpWignerD {
public:
  //============ CONSTRUCTION ==================
  AmpWignerD(ComPWA::Spin spin = ComPWA::Spin(0), unsigned int mu = 0,
             unsigned int muPrime = 0);

  virtual ~AmpWignerD(){};

  static std::shared_ptr<AmpWignerD>
  Factory(std::shared_ptr<PartList> partL,
          const boost::property_tree::ptree &pt);

  //================ EVALUATION =================

  virtual double evaluate(const ComPWA::DataPoint &point, int pos1,
                          int pos2) const;

  static double dynamicalFunction(ComPWA::Spin J, ComPWA::Spin mu,
                                  ComPWA::Spin muPrime, double cosTheta);

  static std::complex<double> dynamicalFunction(double cosAlpha, double cosBeta,
                                                double cosGamma, ComPWA::Spin J,
                                                ComPWA::Spin mu,
                                                ComPWA::Spin muPrime);

  //============ SET/GET =================

  void SetSpin(ComPWA::Spin s) { _spin = s; }
  ComPWA::Spin GetSpin() const { return _spin; }

  void SetMu(ComPWA::Spin s) { _mu = s; }
  ComPWA::Spin GetMu() const { return _mu; }

  void SetMuPrime(ComPWA::Spin s) {
    _helicities = std::pair<ComPWA::Spin, ComPWA::Spin>(0, s);
  }
  ComPWA::Spin GetMuPrime() const {
    return _helicities.first - _helicities.second;
  }

  void SetHelicities(std::pair<ComPWA::Spin, ComPWA::Spin> hel) {
    _helicities = hel;
  }

  std::pair<ComPWA::Spin, ComPWA::Spin> GetHelicities() const {
    return _helicities;
  }

  //=========== FUNCTIONTREE =================

  virtual std::shared_ptr<ComPWA::FunctionTree>
  GetTree(const ComPWA::ParameterList &sample, int posTheta, int posPhi,
          std::string suffix = "");

protected:
  ComPWA::Spin _spin;
  ComPWA::Spin _mu;
  std::pair<ComPWA::Spin, ComPWA::Spin> _helicities;
};

class WignerDStrategy : public Strategy {
public:
  WignerDStrategy(const std::string resonanceName)
      : Strategy(ParType::MDOUBLE), name(resonanceName) {}

  virtual const std::string to_str() const { return ("WignerD of " + name); }

  virtual void execute(ParameterList &paras,
                       std::shared_ptr<Parameter> &out);

protected:
  std::string name;
};

} // namespace AmplitudeSum
} // namespace Physics
} // namespace ComPWA

#endif
