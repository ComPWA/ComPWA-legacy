// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "../RootIO/RootEfficiency.hpp"

#include "TH2.h"

#include "Core/Event.hpp"
#include "Core/Exceptions.hpp"

namespace ComPWA {
namespace Data {

RootEfficiency::RootEfficiency(TEfficiency *eff)
    : effHist(new TEfficiency(*eff)) {
  LOG(DEBUG) << "RootEfficiency: creating efficiency from existing "
                "TEfficiency object!";
}

RootEfficiency::RootEfficiency(TH1 *passed, TH1 *total)
    : effHist(new TEfficiency(*passed, *total)) {
  LOG(DEBUG) << "RootEfficiency: creating efficiency from two TH2D objects!";
}

RootEfficiency::RootEfficiency(const RootEfficiency &) {}

double RootEfficiency::evaluate(const DataPoint &point) const {
  //	double m13sq = point.getVal("m13sq");
  //	double m23sq = point.getVal("m23sq");
  double m13sq = point.KinematicVariableList[1];
  double m23sq = point.KinematicVariableList[0];

  TH2D *test = (TH2D *)effHist->GetPassedHistogram();
  int globalBin = test->FindBin(m23sq, m13sq);
  return effHist->GetEfficiency(globalBin);
}

// ------------------ ROOTANGLEEFFICIENCY ---------------------

double RootAngleEfficiency::evaluate(const DataPoint &point) const {
  double m23sq = point.KinematicVariableList[0];
  double angle = point.KinematicVariableList[8];

  TH2D *test = (TH2D *)effHist->GetPassedHistogram();
  int globalBin = test->FindBin(m23sq, angle);
  return effHist->GetEfficiency(globalBin);
}

} // namespace Data
} /* namespace ComPWA */
