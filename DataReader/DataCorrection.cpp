// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/Kinematics.hpp"
#include "DataReader/DataCorrection.hpp"

namespace ComPWA {

MomentumCorrection::MomentumCorrection(
    std::vector<ComPWA::DataReader::CorrectionTable> inCorr, std::string t)
    : Corrections(inCorr), Title(t) {}

double MomentumCorrection::correction(Event &ev) {
  double w = 1;
  for (int i = 0; i < ev.numParticles(); i++) {
    Particle p = ev.particle(i);
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
  LOG(info) << "MomentumCorrection::Print() | " << Title;
  for (int i = 0; i < Corrections.size(); i++)
    Corrections.at(i).Print();
  return;
}
}
