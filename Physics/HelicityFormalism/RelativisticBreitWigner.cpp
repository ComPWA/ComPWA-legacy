//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//     Peter Weidenkaff - correct nominator, using dataPoint for data handling
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>

#include "boost/property_tree/ptree.hpp"

#include "Core/DataPointStorage.hpp"
#include "Core/Kinematics.hpp"
#include "Physics/HelicityFormalism/RelativisticBreitWigner.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

std::complex<double> RelativisticBreitWigner::Evaluate(const dataPoint &point,
                                                       int pos) const {
  std::complex<double> result = dynamicalFunction(
      point.getVal(pos), _mass->GetValue(), _massA, _massB, _width->GetValue(),
      (double)_spin, _mesonRadius->GetValue(), _ffType);
  return result;
}

std::complex<double> RelativisticBreitWigner::dynamicalFunction(
    double mSq, double mR, double ma, double mb, double width, unsigned int J,
    double mesonRadius, formFactorType ffType) {

  std::complex<double> i(0, 1);
  double sqrtS = sqrt(mSq);

  double barrier =
      Kinematics::FormFactor(sqrtS, ma, mb, J, mesonRadius, ffType) /
      Kinematics::FormFactor(mR, ma, mb, J, mesonRadius, ffType);

  std::complex<double> qTerm = std::pow((Kinematics::phspFactor(sqrtS, ma, mb) /
                                         Kinematics::phspFactor(mR, ma, mb)) *
                                            mR / sqrtS,
                                        (2 * J + 1));

  // Calculate coupling constant to final state
  //  std::complex<double> g_final =
  //      widthToCoupling(mSq, mR, width, ma, mb, J, mesonRadius, ffType);
  std::complex<double> g_final = 1.0;

  /*Coupling constant from production reaction. In case of a particle decay
   * the production coupling doesn't depend in energy since the CM energy
   * is in the (RC) system fixed to the mass of the decaying particle */
  double g_production = 1;

  std::complex<double> denom = std::complex<double>(mR * mR - mSq, 0) +
                               (-1.0) * i * sqrtS * (width * qTerm * barrier);

  std::complex<double> result = g_final * g_production / denom;

#ifndef NDEBUG
  if (std::isnan(result.real()) || std::isnan(result.imag())) {
    std::cout << "RelativisticBreitWigner::dynamicalFunction() | " << barrier
              << " " << mR << " " << mSq << " " << ma << " " << mb << std::endl;
    return 0;
  }
#endif

  return result;
}

std::shared_ptr<RelativisticBreitWigner>
RelativisticBreitWigner::Factory(const boost::property_tree::ptree &pt) {
  LOG(trace) << "RelativisticBreitWigner::Factory() | Construction....";
  auto obj = std::make_shared<RelativisticBreitWigner>();

  int id = pt.get<double>("DecayParticle.<xmlattr>.Id");
  auto partProp = PhysConst::Instance()->FindParticle(id);
  obj->SetMass(std::make_shared<DoubleParameter>(partProp.GetMassPar()));

  auto decayTr = partProp.GetDecayInfo();
  if (partProp.GetDecayType() != "relativisticBreitWigner")
    throw std::runtime_error(
        "RelativisticBreitWigner::Factory() | Decay type does not match! ");

  auto width = ComPWA::DoubleParameterFactory(decayTr.get_child("Width"));
  obj->SetWidth(std::make_shared<DoubleParameter>(width));
  auto mesonRadius =
      ComPWA::DoubleParameterFactory(decayTr.get_child("MesonRadius"));
  obj->SetMesonRadius(std::make_shared<DoubleParameter>(mesonRadius));

  // Get masses of decay products
  auto decayProducts = pt.get_child("DecayProducts");
  std::vector<int> daughterId;
  for (auto i : decayProducts) {
    daughterId.push_back(i.second.get<int>("<xmlattr>.Id"));
  }
  if (daughterId.size() != 2)
    throw boost::property_tree::ptree_error(
        "AmpWignerD::Factory() | Expect exactly two decay products (" +
        std::to_string(decayProducts.size()) + " given)!");

  obj->SetDecayMassA(
      PhysConst::Instance()->FindParticle(daughterId.at(0)).GetMass());
  obj->SetDecayMassB(
      PhysConst::Instance()->FindParticle(daughterId.at(1)).GetMass());

  LOG(trace)
      << "RelativisticBreitWigner::Factory() | Construction of the decay "
      << partProp.GetName() << " -> "
      << PhysConst::Instance()->FindParticle(daughterId.at(0)).GetName()
      << " + "
      << PhysConst::Instance()->FindParticle(daughterId.at(1)).GetName();

  return obj;
}

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */
