// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Data/DataCorrection.hpp"

#include "Core/Kinematics.hpp"

namespace ComPWA {
namespace Data {

MomentumCorrection::MomentumCorrection(
    std::vector<ComPWA::Data::CorrectionTable> inCorr, std::string t)
    : Corrections(inCorr), Title(t) {}

double MomentumCorrection::correction(Event &ev) {
  double w = 1;
  for (unsigned int i = 0; i < ev.ParticleList.size(); ++i) {
    Particle p = ev.ParticleList[i];
    int charge = p.charge();
    double mom = p.threeMomentum();
    double corr;
    try {
      corr = Corrections.at(i).GetValue(charge, mom) + 1;
    } catch (std::exception &ex) {
      throw std::runtime_error(
          "MomentumCorrection::getCorrection() | Number of histograms is "
          "expected to be number of final state particles!");
    }
    w *= corr;
  }
  return w;
}

void MomentumCorrection::print() const {
  LOG(INFO) << "MomentumCorrection::Print() | " << Title;
  for (unsigned int i = 0; i < Corrections.size(); ++i)
    Corrections.at(i).Print();
  return;
}
} // namespace Data
} // namespace ComPWA
