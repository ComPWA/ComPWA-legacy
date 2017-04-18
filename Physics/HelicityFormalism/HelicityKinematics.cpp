//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

#include "Core/Event.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Particle.hpp"
#include "Core/PhysConst.hpp"

#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/qft++/Vector4.h"
#include "Tools/RootGenerator.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

HelicityKinematics::HelicityKinematics(std::vector<int> initialState,
                                       std::vector<int> finalState)
    : Kinematics(initialState, finalState) {
  assert(_initialState.size() == 1);
  _idMother = _initialState.at(0);
  auto motherProp = PhysConst::Instance()->FindParticle(_idMother);
  _nameMother = motherProp.GetName();
  _M = motherProp.GetMass();
  _Msq = _M * _M;
  _spinM = motherProp.GetSpin();

  // Creating unique title
  std::stringstream stream;
  stream << "(";
  for (auto i : _initialState)
    stream << std::to_string(i) << " ";
  stream << ")->(";
  for (auto i : _finalState)
    stream << std::to_string(i) << " ";
  stream << ")";

  LOG(info) << "HelicityKinematics::HelicityKinematics() | Initialize reaction "
            << stream.str();
}

HelicityKinematics::HelicityKinematics(boost::property_tree::ptree pt) {

  auto initialS = pt.get_child("HelicityKinematics.InitialState");
  for (auto i : initialS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = PhysConst::Instance()->FindParticle(name);
    _initialState.push_back(partP.GetId());
  }

  auto finalS = pt.get_child("HelicityKinematics.FinalState");
  for (auto i : finalS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = PhysConst::Instance()->FindParticle(name);
    _finalState.push_back(partP.GetId());
  }

  // Creating unique title
  std::stringstream stream;
  stream << "(";
  for (auto i : _initialState)
    stream << std::to_string(i) << " ";
  stream << ")->(";
  for (auto i : _finalState)
    stream << std::to_string(i) << " ";
  stream << ")";

  LOG(info) << "HelicityKinematics::HelicityKinematics() | Initialize reaction "
            << stream.str();

  assert(_initialState.size() == 1);
  _idMother = _initialState.at(0);
  auto motherProp = PhysConst::Instance()->FindParticle(_idMother);
  _nameMother = motherProp.GetName();
  _M = motherProp.GetMass();
  _Msq = _M * _M;
  _spinM = motherProp.GetSpin();
}

HelicityKinematics::~HelicityKinematics() {}

bool HelicityKinematics::IsWithinPhsp(const dataPoint &point) const {
  // TODO: implementation

  return true;
}

std::pair<double, double> HelicityKinematics::GetInvMassBounds(SubSystem sys) {

  /* We use the formulae from (PDG2016 Kinematics Fig.47.3). I hope the
   * generalization to n-body decays is correct.
   */
  double min = 0; // min = (m1+m2)^2 (Dalitz decay)
  for (auto i : sys.GetFinalStateA())
    min += PhysConst::Instance()->FindParticle(_finalState.at(i)).GetMass();
  for (auto i : sys.GetFinalStateB())
    min += PhysConst::Instance()->FindParticle(_finalState.at(i)).GetMass();
  min = min * min;

  double max = _M; // max = (M-m3)^2 (Dalitz decay)
  for (auto i : sys.GetRecoilState())
    max -= PhysConst::Instance()->FindParticle(_finalState.at(i)).GetMass();
  max = max * max;

  LOG(trace) << "HelicityKinematics::GetInvMassBounds() | Bounds for SubSystem "
             << sys << " are [" << min << "," << max << "]";
  return std::pair<double, double>(min, max);
}

double HelicityKinematics::calculatePSArea() {
  double result(0);
  double precision(1); // relative uncertainty
  double weights(0);
  double precisionLimit(1e-5);
  unsigned int maxCalls(1e6);
  unsigned int call(0);

  std::shared_ptr<ComPWA::Generator> gen(new ComPWA::Tools::RootGenerator(0));

  while (precision > precisionLimit && call < maxCalls) {
    precision = result;
    ComPWA::Event ev;
    gen->generate(ev);
    weights += ev.getWeight();
    call++;
    result = weights / call;
    precision = std::fabs(result - precision) / result;
  }

  LOG(debug) << "HelicityKinematics::calculatePSarea() | "
             << "(" << result << "+-" << precision * result
             << ") GeV^4 relAcc [%]: " << 100 * precision;

  return result;
}

void HelicityKinematics::EventToDataPoint(const Event &event,
                                          dataPoint &point) const {

  for (auto i : _listSubSystem) {
    QFT::Vector4<double> recoil, finalA, finalB;
    for (auto s : i.GetRecoilState()) {
      recoil += QFT::Vector4<double>(event.getParticle(s).getFourMomentum());
    }
    for (auto s : i.GetFinalStateA()) {
      finalA += QFT::Vector4<double>(event.getParticle(s).getFourMomentum());
    }
    for (auto s : i.GetFinalStateB()) {
      finalB += QFT::Vector4<double>(event.getParticle(s).getFourMomentum());
    }
    QFT::Vector4<double> dd = finalA + finalB;
    double mSq = dd.Mass2();

    finalA.Boost(dd);
    finalA.Rotate(dd.Phi(), dd.Theta(), (-1) * dd.Phi());
    double cosTheta = finalA.CosTheta();
    double phi = finalA.Phi();
    if (cosTheta > 1 || cosTheta < -1 || phi > M_PI || phi < (-1) * M_PI ||
        std::isnan(cosTheta) || std::isnan(phi)) {
      throw BeyondPhsp("HelicityKinematics::EventToDataPoint() |"
                       " Point beypond phase space boundaries!");
    }

    // LOG(trace)<<dd.Mass2()<< " " <<std::cos(finalA.Theta())<<"
    // "<<finalA.Phi();

    point.GetPoint().push_back(mSq);
    point.GetPoint().push_back(cosTheta);
    point.GetPoint().push_back(phi);
  }
}

double HelicityKinematics::FormFactor(double sqrtS, double ma, double mb,
                                      double spin, double mesonRadius,
                                      formFactorType type) {
  if (type == formFactorType::BlattWeisskopf && spin == 0) {
    return 1.0;
  }

  std::complex<double> qValue = Kinematics::qValue(sqrtS, ma, mb);
  return HelicityKinematics::FormFactor(sqrtS, ma, mb, spin, mesonRadius,
                                        qValue, type);
}

double HelicityKinematics::FormFactor(double sqrtS, double ma, double mb,
                                      double spin, double mesonRadius,
                                      std::complex<double> qValue,
                                      formFactorType type) {
  if (mesonRadius == 0)
    return 1; // disable form factors
  if (type == formFactorType::noFormFactor)
    return 1; // disable form factors
  if (type == formFactorType::BlattWeisskopf && spin == 0) {
    return 1.0;
  }

  // From factor for a0(980) used by Crystal Barrel Phys.Rev.D78-074023
  if (type == formFactorType::CrystalBarrel) {
    if (spin == 0) {
      double qSq = std::norm(qValue);
      double alpha = mesonRadius * mesonRadius / 6;
      return std::exp(-alpha * qSq);
    } else
      throw std::runtime_error("Kinematics::FormFactor() | "
                               "Form factors of type " +
                               std::string(formFactorTypeString[type]) +
                               " are implemented for spin 0 only!");
  }

  // Blatt-Weisskopt form factors with normalization F(x=mR) = 1.
  // Reference: S.U.Chung Annalen der Physik 4(1995) 404-430
  // z = q / (interaction range). For the interaction range we assume
  // 1/mesonRadius
  if (type == formFactorType::BlattWeisskopf) {
    if (spin == 0)
      return 1;
    double qSq = std::norm(qValue);
    double z = qSq * mesonRadius * mesonRadius;
    /* Events below threshold
     * What should we do if event is below threshold? Shouldn't really influence
     * the result
     * because resonances at threshold don't have spin(?) */
    z = std::fabs(z);

    if (spin == 1) {
      return (sqrt(2 * z / (z + 1)));
    } else if (spin == 2) {
      return (sqrt(13 * z * z / ((z - 3) * (z - 3) + 9 * z)));
    } else if (spin == 3) {
      return (
          sqrt(277 * z * z * z / (z * (z - 15) * (z - 15) + 9 * (2 * z - 5))));
    } else if (spin == 4) {
      return (sqrt(12746 * z * z * z * z /
                   ((z * z - 45 * z + 105) * (z * z - 45 * z + 105) +
                    25 * z * (2 * z - 21) * (2 * z - 21))));
    } else
      throw std::runtime_error(
          "Kinematics::FormFactor() | Form factors of type " +
          std::string(formFactorTypeString[type]) +
          " are implemented for spins up to 4!");
  }
  throw std::runtime_error("Kinematics::FormFactor() | Form factor type " +
                           std::to_string((long long int)type) +
                           " not specified!");

  return 0;
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
