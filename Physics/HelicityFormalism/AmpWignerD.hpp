// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Contains AmpWignerD class which provides WignerD functions.
///

#ifndef COMPWA_PHYSICS_HELICITY_FORMALISM_AMPWIGNERD_HPP_
#define COMPWA_PHYSICS_HELICITY_FORMALISM_AMPWIGNERD_HPP_

#include <memory>
#include <vector>

#include "Core/FunctionTree.hpp"
#include "Core/Functions.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Properties.hpp"
#include "Core/Spin.hpp"

namespace ComPWA {

struct DataPoint;

namespace Physics {
namespace HelicityFormalism {

///
/// \class AmpWignerD
/// Angular distribution based on WignerD functions
///
class AmpWignerD {
public:
  //============ CONSTRUCTION ==================
  AmpWignerD(ComPWA::Spin spin, ComPWA::Spin muPrime, ComPWA::Spin mu);

  virtual ~AmpWignerD(){};
  //================ EVALUATION =================

  virtual std::complex<double> evaluate(const ComPWA::DataPoint &point,
                                        int pos1, int pos2) const;

  static double dynamicalFunction(ComPWA::Spin J, ComPWA::Spin muPrime,
                                  ComPWA::Spin mu, double beta);

  static std::complex<double> dynamicalFunction(ComPWA::Spin J,
                                                ComPWA::Spin muPrime,
                                                ComPWA::Spin mu, double alpha,
                                                double beta, double gamma);

  //=========== FUNCTIONTREE =================
  virtual std::shared_ptr<ComPWA::FunctionTree>
  tree(const ComPWA::ParameterList &sample, int posTheta, int posPhi,
       std::string suffix = "");

protected:
  ComPWA::Spin J;
  ComPWA::Spin MuPrime;
  ComPWA::Spin Mu;
};

class WignerDStrategy : public Strategy {
public:
  WignerDStrategy(const std::string resonanceName)
      : Strategy(ParType::MCOMPLEX), name(resonanceName) {}

  virtual const std::string to_str() const { return ("WignerD of " + name); }

  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);

protected:
  std::string name;
};

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA

#endif
