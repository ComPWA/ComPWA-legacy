
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
#include <algorithm>

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
    stream << PhysConst::Instance()->FindParticle(i).GetName() << " ";
  stream << ")->( ";
  for (auto i : _finalState)
    stream << PhysConst::Instance()->FindParticle(i).GetName() << " ";
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
  stream << "( ";
  for (auto i : _initialState)
    stream << PhysConst::Instance()->FindParticle(i).GetName() << " ";
  stream << ")->( ";
  for (auto i : _finalState)
    stream << PhysConst::Instance()->FindParticle(i).GetName() << " ";
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
  int subSystemID = 0;
  int pos = 0;
  while ((pos + 2) < point.Size()) {
    auto invMassBounds = GetInvMassBounds(GetSubSystem(subSystemID));
    if (point.GetValue(pos) < invMassBounds.first ||
        point.GetValue(pos) > invMassBounds.second)
      return false;
    if (point.GetValue(pos + 1) < -1 || point.GetValue(pos + 1) > +1)
      return false;
    if (point.GetValue(pos + 2) < 0 || point.GetValue(pos + 2) > 2 * M_PI)
      return false;

    pos += 3;
    subSystemID++;
  }

  return true;
}

std::pair<double, double>
HelicityKinematics::GetInvMassBounds(SubSystem sys) const {

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

  //  LOG(trace) << "HelicityKinematics::GetInvMassBounds() | Bounds for
  //  SubSystem "
  //             << sys << " are [" << min << "," << max << "]";
  return std::pair<double, double>(min, max);
}

double HelicityKinematics::EventDistance(Event &evA, Event &evB) const {
  dataPoint pointA, pointB;
  std::vector<int> finalAll;
  for (int i = 0; i < _finalState.size(); i++) {
    finalAll.push_back(i);
  }

  /* We fill dataPoints with variables from a set of SubSystems such that
   * each final state momentum is used in the calculation.
   * E.g. for a three particle decay (Dalitz plot) with final state particles
   * (123) we use (12) and (13). The invariant masses are used
   * in the next step to calculate the distance in phase space. We hope that
   * this
   * set of invariant masses specifies the position in phase space uniquely for
   * all reactions.
   */
  for (int i = 1; i < finalAll.size(); i++) {
    std::vector<int> finalA;
    finalA.push_back(0);
    std::vector<int> finalB;
    finalB.push_back(i);
    std::vector<int> recoil(finalAll);
    recoil.erase(std::remove(recoil.begin(), recoil.end(), 0), recoil.end());
    recoil.erase(std::remove(recoil.begin(), recoil.end(), i), recoil.end());

    SubSystem sys(recoil, finalA, finalB);
    EventToDataPoint(evA, pointA, sys);
    EventToDataPoint(evB, pointB, sys);
  }

  /* We assume that variables in dataPoint are orders like (m1,theta1,phi1),
   * (m2,theta2,phi2),... Therefore, the dataPoint size needs to be a multiple
   * of 3. */
  assert(pointA.Size() % 3 == 0);

  double distSq = 0;
  int pos = 0;
  while (pos < pointA.Size()) {
    double dist = (pointA.GetValue(pos) - pointB.GetValue(pos));
    distSq += dist * dist;
    pos += 2;
  }
  return std::sqrt(distSq);
}

/* This is not working since I have no idea how we can check if a dataPoint is
 * within phase space boundaries if we have only the invariant masses.

double HelicityKinematics::calculatePSArea() {
int precision = 2000000; // sample size
auto gen = ComPWA::Tools::RootGenerator(0);

std::vector<int> finalAll;
for (int i = 0; i < _finalState.size(); i++) {
  finalAll.push_back(i);
}
std::vector<std::pair<double, double>> limits;
double totalVolume = 1;

* We fill dataPoints with variables from a set of SubSystems such that
 * each final state momentum is used in the calculation.
 * E.g. for a three particle decay (Dalitz plot) with final state particles
 * (123) we use (12) and (13). The invariant masses are used
 * in the next step to calculate the distance in phase space. We hope that
 * this
 * set of invariant masses specifies the position in phase space uniquely for
 * all reactions.
 *
for (int i = 1; i < finalAll.size(); i++) {
  std::vector<int> finalA;
  finalA.push_back(0);
  std::vector<int> finalB;
  finalB.push_back(i);
  std::vector<int> recoil(finalAll);
  recoil.erase(std::remove(recoil.begin(), recoil.end(), 0), recoil.end());
  recoil.erase(std::remove(recoil.begin(), recoil.end(), i), recoil.end());

  SubSystem sys(recoil, finalA, finalB);
  GetDataID(sys);
  auto bounds = GetInvMassBounds(sys);
  limits.push_back(bounds);
  totalVolume *= (bounds.second - bounds.first);
}
int inside = 0;
for (int p = 0; p < precision; p++) {
  dataPoint point;
  for (int i = 0; i < limits.size(); i++) {
    double mSq = gen.GetUniform(limits.at(i).first, limits.at(i).second);
    mSq = gen.GetUniform(1,2);
    point.GetPoint().push_back(mSq); // mass
    point.GetPoint().push_back(0);   // cosTheta
    point.GetPoint().push_back(1.0); // phi
  }

  if (IsWithinPhsp(point))
    inside++;
}
std::cout<<totalVolume<<" "<<precision<<" "<<inside<<std::endl;
double phspVolume = inside / (double)precision * totalVolume;
LOG(info) << "HelicityKinematics::calculatePSArea() | Phase space volume: "
<< phspVolume;
return phspVolume;
}
*/

double HelicityKinematics::calculatePSArea() {
  int precision = 500; // sample size
  int numNearestN = 10; //Number of nearest neighbours

  // Generate phase space sample
  auto gen = ComPWA::Tools::RootGenerator(123456);
  auto sample = ComPWA::DataReader::RootReader();
  int m = 0;
  while (m < precision) {
    Event tmp;
    gen.Generate(tmp);
    double ampRnd = gen.GetUniform(0,1);
    if (ampRnd > tmp.getWeight())
      continue;
    m++;
    sample.pushEvent(tmp);
  }

  std::vector<double> localDistance(sample.getNEvents(),0.0);
  for (int n = 0; n < sample.getNEvents(); n++) {
    std::vector<double> dist;
    for (int m = 0; m < sample.getNEvents(); m++) {
      double d = EventDistance(sample.getEvent(n), sample.getEvent(m));
      dist.push_back(d);
    }
    std::sort(dist.begin(), dist.end()); // ascending
    for(int i = 0; i<localDistance.size(); i++)
      localDistance.at(i) += dist.at(i);
  }
  for(int i=0; i<localDistance.size(); i++)
    localDistance.at(i) /= sample.getNEvents();

  int nDim = _finalState.size() - 1;
  
  std::vector<double> avgVecVol(sample.getNEvents());
  for(int i=0; i<localDistance.size(); i++)
    avgVecVol.at(i) = std::pow(M_PI, nDim / 2) / std::tgamma(nDim / 2 + 1) *
                     std::pow(localDistance.at(i), nDim);
  
  std::vector<double> avgVecDensity(sample.getNEvents());
  for(int i=0; i<localDistance.size(); i++)
    avgVecDensity.at(i) = sample.getNEvents() / (i / avgVecVol.at(i));
  
  std::cout<<"Average distance | volume | density:"<<std::endl;
  for(int i=0; i<precision/4; i++)
    std::cout<<i<<" "<<localDistance.at(i)<<" "<<avgVecVol.at(i)<<" "<<avgVecDensity.at(i)<<std::endl;
  
  double phspVolume = avgVecDensity.at(numNearestN);
  LOG(info) << "HelicityKinematics::calculatePSArea() | Phase space volume: "
            << phspVolume;
  return phspVolume;
}

void HelicityKinematics::EventToDataPoint(const Event &event,
                                          dataPoint &point) const {
  for (auto i : _listSubSystem)
    EventToDataPoint(event, point, i);
  return;
}

void HelicityKinematics::EventToDataPoint(const Event &event, dataPoint &point,
                                          SubSystem sys) const {
  QFT::Vector4<double> recoil, finalA, finalB;
  for (auto s : sys.GetRecoilState()) {
    recoil += QFT::Vector4<double>(event.getParticle(s).getFourMomentum());
  }
  for (auto s : sys.GetFinalStateA()) {
    finalA += QFT::Vector4<double>(event.getParticle(s).getFourMomentum());
  }
  for (auto s : sys.GetFinalStateB()) {
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
