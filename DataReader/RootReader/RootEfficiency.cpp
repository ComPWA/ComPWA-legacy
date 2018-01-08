// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/Exceptions.hpp"
#include "Core/DataPoint.hpp"
#include "DataReader/RootReader/RootEfficiency.hpp"

namespace ComPWA {
namespace DataReader {

RootEfficiency::RootEfficiency(TEfficiency *eff)
    : effHist(new TEfficiency(*eff)) {
  LOG(debug) << "RootEfficiency: creating efficiency from existing "
                "TEfficiency object!";
}
  
RootEfficiency::RootEfficiency(TH1 *passed, TH1 *total)
    : effHist(new TEfficiency(*passed, *total)) {
  LOG(debug)
      << "RootEfficiency: creating efficiency from two TH2D objects!";
}
  
RootEfficiency::RootEfficiency(const RootEfficiency &) {}
  
double RootEfficiency::evaluate(const DataPoint &point) const {
  //	double m13sq = point.getVal("m13sq");
  //	double m23sq = point.getVal("m23sq");
  double m13sq = point.value(1);
  double m23sq = point.value(0);

  TH2D *test = (TH2D *)effHist->GetPassedHistogram();
  int globalBin = test->FindBin(m23sq, m13sq);
  return effHist->GetEfficiency(globalBin);
}
  
  
  // ------------------ ROOTANGLEEFFICIENCY ---------------------

double RootAngleEfficiency::evaluate(const DataPoint &point) const {
  double m23sq = point.value(0);
  double angle = point.value(8);

  TH2D *test = (TH2D *)effHist->GetPassedHistogram();
  int globalBin = test->FindBin(m23sq, angle);
  return effHist->GetEfficiency(globalBin);
}

} /* namespace DataReader */
} /* namespace ComPWA */
