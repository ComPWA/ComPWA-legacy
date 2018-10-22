// Copyright (c) 2013 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

#include "Estimator/ChiOneD/ChiOneD.hpp"

#include "Data/Data.hpp"
#include "Core/ParameterList.hpp"
#include "Core/AmpIntensity.hpp"
#include "Core/Event.hpp"

namespace ComPWA {

namespace Estimator {
namespace ChiOneD {

ChiOneD::ChiOneD(std::shared_ptr<Kinematics> kin,
                 std::shared_ptr<AmpIntensity> inPIF,
                 std::shared_ptr<ComPWA::Data::Data> inDIF)
    : pPIF_(inPIF), pDIF_(inDIF) {}

double ChiOneD::controlParameter(ParameterList &minPar) {
  unsigned int nBins = pDIF_->numBins();

  double chi = 0;
  for (unsigned int bin = 0; bin < nBins; bin++) {
    double m12, weight;
    pDIF_->bin(bin, m12, weight);

    std::vector<double> x;
    x.push_back(m12);
    double intens;
    //TODO: How do we solve this?
    //intens = pPIF_->Intensity(x);
    assert( 0 && "ChiOneD::controlParameter | Not implemented!");

    chi += (weight - intens) * (weight - intens) /
           (double)nBins; // Just for binned data
  }

  return chi;
}

} /* namespace ChiOneD */
} /* namespace Estimator */
} /* namespace ComPWA */
