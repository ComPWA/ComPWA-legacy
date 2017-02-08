 
            /*
 * DataCorrection.cpp
 *
 *  Created on: Aug 27, 2015
 *      Author: weidenka
 */

#include "Core/Kinematics.hpp"
#include "DataReader/DataCorrection.hpp"

namespace ComPWA {

MomentumCorrection::MomentumCorrection(std::vector<CorrectionTable> inCorr,
                                       std::string t)
    : corrections(inCorr), title(t) {
  if (corrections.size() != Kinematics::instance()->GetNumberOfParticles())
    throw std::runtime_error(
        "RootCorrection::RootCorrection() | Number of histograms is "
        "expected to be number of final state particles!");
}

double MomentumCorrection::getCorrection(Event &ev) {
  double w = 1;
  for (int i = 0; i < ev.getNParticles(); i++) {
    Particle p = ev.getParticle(i);
    int charge = p.getCharge();
    double mom = p.getThreeMomentum();
    double corr;
    try {
      corr = corrections.at(i).GetValue(charge, mom) + 1;
    } catch (...) { // if no correction value is available we set it to one
      corr = 1.0;
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
