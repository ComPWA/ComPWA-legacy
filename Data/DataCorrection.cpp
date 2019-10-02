// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Data/DataCorrection.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Properties.hpp"

namespace ComPWA {
namespace Data {

MomentumCorrection::MomentumCorrection(
    ComPWA::ParticleList PartList_,
    std::vector<ComPWA::Data::CorrectionTable> inCorr, std::string t)
    : PartList(PartList_), Corrections(inCorr), Title(t) {}

double MomentumCorrection::correction(Event &ev) const {
  double w = 1;
  for (unsigned int i = 0; i < ev.ParticleList.size(); ++i) {
    Particle p = ev.ParticleList[i];
    int charge =
        findParticle(PartList, p.pid()).getQuantumNumber<int>("charge");
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

} // namespace Data
} // namespace ComPWA
