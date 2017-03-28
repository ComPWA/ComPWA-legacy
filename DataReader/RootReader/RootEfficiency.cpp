//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff -
//-------------------------------------------------------------------------------
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
  
double RootEfficiency::Evaluate(const dataPoint &point) const {
  //	double m13sq = point.getVal("m13sq");
  //	double m23sq = point.getVal("m23sq");
  double m13sq = point.GetValue(1);
  double m23sq = point.GetValue(0);

  TH2D *test = (TH2D *)effHist->GetPassedHistogram();
  int globalBin = test->FindBin(m23sq, m13sq);
  return effHist->GetEfficiency(globalBin);
}
  
  
  // ------------------ ROOTANGLEEFFICIENCY ---------------------

double RootAngleEfficiency::Evaluate(const dataPoint &point) const {
  double m23sq = point.GetValue(0);
  double angle = point.GetValue(8);

  TH2D *test = (TH2D *)effHist->GetPassedHistogram();
  int globalBin = test->FindBin(m23sq, angle);
  return effHist->GetEfficiency(globalBin);
}

} /* namespace DataReader */
} /* namespace ComPWA */
