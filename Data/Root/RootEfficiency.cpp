// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Data/Root/RootEfficiency.hpp"
#include "Core/Event.hpp"
#include "Core/Exceptions.hpp"
#include "Data/DataSet.hpp"

#include "TH2.h"

namespace ComPWA {
namespace Data {
namespace Root {

RootEfficiency::RootEfficiency(TEfficiency *eff)
    : effHist(new TEfficiency(*eff)) {
  LOG(DEBUG) << "RootEfficiency: creating efficiency from existing "
                "TEfficiency object!";
}

RootEfficiency::RootEfficiency(TH1 *passed, TH1 *total)
    : effHist(new TEfficiency(*passed, *total)) {
  LOG(DEBUG) << "RootEfficiency: creating efficiency from two TH2D objects!";
}

std::vector<double> RootEfficiency::evaluate(const DataSet &dataset) const {
  throw std::runtime_error(
      "RootEfficiency::evaluate(): is currently not implemented!");
  return {};
}

// ------------------ ROOTANGLEEFFICIENCY ---------------------

std::vector<double>
RootAngleEfficiency::evaluate(const DataSet &dataset) const {
  throw std::runtime_error(
      "RootAngleEfficiency::evaluate(): is currently not implemented!");
  return {};
}

} // namespace Root
} // namespace Data
} // namespace ComPWA
