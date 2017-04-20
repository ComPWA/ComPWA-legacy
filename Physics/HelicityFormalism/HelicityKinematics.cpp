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

#include <numeric>
#include <cmath>
#include "Core/Event.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Particle.hpp"
#include "Core/PhysConst.hpp"
#include "DataReader/RootReader/RootReader.hpp"

#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/qft++/Vector4.h"
#include "Tools/RootGenerator.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

HelicityKinematics::HelicityKinematics(std::vector<pid> initialState,
                                       std::vector<pid> finalState)
    : Kinematics(initialState, finalState) {
  assert(_initialState.size() == 1);
  pid idMother = _initialState.at(0);
  auto motherProp = PhysConst::Instance()->FindParticle(idMother);
  _M = motherProp.GetMass();
  _Msq = _M * _M;
  _spinM = motherProp.GetSpin();

  // Creating unique title
  std::stringstream stream;
  stream << "( ";
  for (auto i : _initialState)
    stream << std::to_string(i) << " ";
  stream << ")->( ";
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
  pid idMother = _initialState.at(0);
  auto motherProp = PhysConst::Instance()->FindParticle(idMother);
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

double HelicityKinematics::EventDistance(Event &evA, Event &evB) const {
  dataPoint pointA, pointB;
  EventToDataPoint(evA, pointA);
  EventToDataPoint(evB, pointB);

  /* We assume that variables in dataPoint are orders like (m1,theta1,phi1),
   * (m2,theta2,phi2),... Therefore, the dataPoint size needs to be a multiple 
   * of 3. */
  assert( pointA.Size()%3 == 0 );
  
  double distSq = 0;
  int pos = 0;
  while( pos<pointA.Size() ){
    double dist = (pointA.GetValue(pos)-pointB.GetValue(pos));
    distSq += dist*dist;
    pos += 2;
  }
  return std::sqrt(distSq);
}

double HelicityKinematics::calculatePSArea() {
  /* We estimate the volume of the phase space for arbirary decays.
   */

  int precision = 100; // sample size
  int localDensityNumber =
      precision / 10; // number of points within a n-dim sphere

  // Generate phase space sample
  auto gen = ComPWA::Tools::RootGenerator(0);
  auto sample = ComPWA::DataReader::RootReader();
  int i = 0;
  while (i < precision) {
    Event tmp;
    gen.generate(tmp);
    double ampRnd = gen.getUniform();
    if (ampRnd > tmp.getWeight())
      continue;
    i++;
    sample.pushEvent(tmp);
  }

  std::vector<double> localDistance;
  for (int n = 0; n < sample.getNEvents(); n++) {
    std::vector<double> dist;
    for (int m = 0; m < sample.getNEvents(); m++) {
      double d = EventDistance(sample.getEvent(n), sample.getEvent(m));
      dist.push_back(d);
    }
    std::sort(dist.begin(), dist.end()); //ascending
    localDistance.push_back(dist.at(localDensityNumber));
  }

  // Calculate average distance that contains @localDensityNumber events
  double avgDistance =
      std::accumulate(localDistance.begin(), localDistance.end(), 0.0) /
      localDistance.size();

  int nDim = GetIrreducibleSetOfVariables().size();
  // N-dimensional sphere with avgDististance as radius
  double avgVolume = std::pow(M_PI, nDim / 2) / std::tgamma(nDim / 2 + 1) *
                     std::pow(avgDistance, nDim);
  double avgDensity = localDensityNumber / avgVolume;

  double phspVolume = sample.getNEvents() / avgDensity;
  LOG(info) << "HelicityKinematics::calculatePSArea() | Phase space volume: "
            << phspVolume;
  return phspVolume;
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
