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
  AmpWignerD(ComPWA::Spin spin, ComPWA::Spin mu, ComPWA::Spin muPrime);

  virtual ~AmpWignerD(){};
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

  void setSpin(ComPWA::Spin s) { J = s; }
  ComPWA::Spin spin() const { return J; }

  void setMu(ComPWA::Spin s) { Mu = s; }
  ComPWA::Spin mu() const { return Mu; }

  void setMuPrime(ComPWA::Spin s) { MuPrime = s; }

  ComPWA::Spin muPrime() const { return MuPrime; }

  //=========== FUNCTIONTREE =================
  virtual std::shared_ptr<ComPWA::FunctionTree>
  tree(const ComPWA::ParameterList &sample, int posTheta, int posPhi,
          std::string suffix = "");

protected:
  ComPWA::Spin J;
  ComPWA::Spin Mu;
  ComPWA::Spin MuPrime;
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
