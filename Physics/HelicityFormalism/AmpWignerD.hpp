//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel
//     Peter Weidenkaff
//-------------------------------------------------------------------------------

//! Angular distribution based on WignerD functions
/*!
 * @file AmpWignerD.hpp
 *\class AmpWignerD
 *The helicity angle for sub system \_subSys is calculated and the value of the
 *WignerD function is returned
 */

#ifndef AMPWIGNER_D
#define AMPWIGNER_D

#include <vector>
#include <memory>

#include "Core/ParameterList.hpp"
#include "Core/Functions.hpp"
#include "Core/DataPoint.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Spin.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

class AmpWignerD {
public:
  AmpWignerD(ComPWA::Spin spin = ComPWA::Spin(0),
             unsigned int mu = 0, unsigned int muPrime = 0);

  virtual ~AmpWignerD(){};

  virtual double Evaluate(dataPoint &point, int pos1, int pos2);

  static double dynamicalFunction(ComPWA::Spin J, ComPWA::Spin mu,
                                  ComPWA::Spin muPrime, double cosTheta);

  static std::complex<double> dynamicalFunction(double cosAlpha, double cosBeta,
                                                double cosGamma, ComPWA::Spin J,
                                                ComPWA::Spin mu,
                                                ComPWA::Spin muPrime);

  virtual std::shared_ptr<FunctionTree> SetupTree(ParameterList &sample,
                                                  std::string suffix = "");

protected:
  unsigned int _varId;
  ComPWA::Spin _spin;
  unsigned int _mu;
  unsigned int _muPrime;
};

class WignerDStrategy : public Strategy {
public:
  WignerDStrategy(const std::string resonanceName)
      : Strategy(ParType::MDOUBLE), name(resonanceName) {}

  virtual const std::string to_str() const { return ("WignerD of " + name); }

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out);

protected:
  std::string name;
};

} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */

#endif
