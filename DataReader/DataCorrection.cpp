// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/Kinematics.hpp"
#include "DataReader/DataCorrection.hpp"

namespace ComPWA {

MomentumCorrection::MomentumCorrection(std::vector<CorrectionTable> inCorr,
                                       std::string t)
    : corrections(inCorr), title(t) {}

double MomentumCorrection::getCorrection(Event &ev) {
  double w = 1;
  for (int i = 0; i < ev.GetNParticles(); i++) {
    Particle p = ev.GetParticle(i);
    int charge = p.GetCharge();
    double mom = p.GetThreeMomentum();
    double corr;
    try {
      corr = corrections.at(i).GetValue(charge, mom) + 1;
    } catch (std::exception &ex) {
      throw std::runtime_error(
          "RootCorrection::RootCorrection() | Number of histograms is "
          "expected to be number of final state particles!");
    }
    w *= corr;
  }
  return w;
}

void MomentumCorrection::Print() const {
  LOG(info) << "MomentumCorrection::Print() | " << title;
  for (int i = 0; i < corrections.size(); i++)
    corrections.at(i).Print();
  return;
}
}
