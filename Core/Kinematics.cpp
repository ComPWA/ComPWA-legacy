//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//		Peter Weidenkaff -
//-------------------------------------------------------------------------------

#include "Core/Kinematics.hpp"
#include "Core/Exceptions.hpp"
#include "Core/PhysConst.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Event.hpp"

namespace ComPWA {

Kinematics *Kinematics::Instance() {
  if (!_inst) {
    throw std::runtime_error("No instance of Kinematics created! "
                             "Create one first!");
  }

  return Kinematics::_inst;
}

Kinematics *Kinematics::_inst = 0;

double Kinematics::qSqValue(double sqrtS, double ma, double mb) {
  double mapb = ma + mb;
  double mamb = ma - mb;
  double xSq = sqrtS * sqrtS;
  double t1 = xSq - mapb * mapb;
  double t2 = xSq - mamb * mamb;
  return (t1 * t2 / (4 * xSq));
}

std::complex<double> Kinematics::qValue(double sqrtS, double ma, double mb) {
  return Kinematics::phspFactor(sqrtS, ma, mb) * 8.0 * M_PI * sqrtS;
}


std::complex<double> Kinematics::phspFactor(double sqrtS, double ma,
                                            double mb) {
  double s = sqrtS * sqrtS;
  std::complex<double> i(0, 1);

  // == Two types of analytic continuation
  // 1) Complex sqrt
  //	std::complex<double> rhoOld;
  //	rhoOld = sqrt(std::complex<double>(qSqValue(sqrtS,ma,mb))) /
  //(8*M_PI*sqrtS); //PDG definition
  //	rhoOld = sqrt(std::complex<double>(qSqValue(sqrtS,ma,mb))) /
  //(0.5*sqrtS); //BaBar definition
  //	return rhoOld; //complex sqrt

  /* 2) Correct analytic continuation
   * proper analytic continuation (PDG 2014 - Resonances (47.2.2))
   * I'm not sure of this is correct for the case of two different masses ma and
   * mb.
   * Furthermore we divide by the factor 16*Pi*Sqrt[s]). This is more or less
   * arbitrary
   * and not mentioned in the reference, but it leads to a good agreement
   * between both
   * approaches.
   */
  std::complex<double> rho;
  double q = std::sqrt(std::fabs(Kinematics::qSqValue(sqrtS, ma, mb) * 4 / s));
  if (s <= 0) { // below 0
    rho = -q / M_PI * std::log(std::fabs((1 + q) / (1 - q)));
  } else if (0 < s && s <= (ma + mb) * (ma + mb)) { // below threshold
    rho = (-2.0 * q / M_PI * atan(1 / q)) / (i * 16.0 * M_PI * sqrtS);
  } else if ((ma + mb) * (ma + mb) < s) { // above threshold
    rho = (-q / M_PI * log(std::fabs((1 + q) / (1 - q))) + i * q) /
          (i * 16.0 * M_PI * sqrtS);
  } else
    throw std::runtime_error("Kinematics::phspFactor| PhspFactor not "
                             "defined at sqrtS = " +
                             std::to_string((long double)sqrtS));
  
#ifndef NDEBUG
  if (rho.real() != rho.real() || rho.imag() != rho.imag()){
    throw std::runtime_error("Kinematics::phspFactor| Result invalid (NaN)!");
  }
#endif

  return rho; // correct analytical continuation
}
  

//! calculated the PHSP volume of the current decay by MC integration
double Kinematics::GetPhspVolume() {
  if (!is_PS_area_calculated_) {
    PS_area_ = calculatePSArea();
    is_PS_area_calculated_ = true;
  }
  return PS_area_;
}

void Kinematics::SetPhspVolume(double vol) {
  PS_area_ = vol;
  is_PS_area_calculated_ = true;
  
  LOG(info)<<"Kinematics::SetPhspVolume() | Setting phase space "
  "volume to "<<std::to_string(GetPhspVolume())<<".";
}

} /* namespace ComPWA */
