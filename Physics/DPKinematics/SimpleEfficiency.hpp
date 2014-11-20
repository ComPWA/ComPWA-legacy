//-------------------------------------------------------------------------------
// Copyright (c) 2014 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff -
//-------------------------------------------------------------------------------
/*! \class SimpleEfficiency
 * @file SimpleEfficiency.hpp
 * Implementation of an efficiency class. Similar functionality to TEfficiency
 */

#ifndef SIMPLEEFFICIENCY_HPP_
#define SIMPLEEFFICIENCY_HPP_
#include "TH1.h"
#include <vector>
#include <memory>

class SimpleEfficiency : public TNamed
{
public:
	~SimpleEfficiency();
	SimpleEfficiency(const TH1& passed, const TH1& total);
	SimpleEfficiency(TString name, TString title, const TH1& passed, const TH1& total);
	SimpleEfficiency():passedHist(NULL),totalHist(NULL),effHist(NULL){};
	SimpleEfficiency(const SimpleEfficiency&);
	double GetTotalEfficiency() { return totalEff; };
	double GetTotalEfficiencyError() { return totalEff_error; };
	double GetEfficiency(int globalBin);
	double GetEfficiencyError(int globalBin);
	double GetEfficiencyErrorLow(int globalBin) { return GetEfficiencyError(globalBin); };
	double GetEfficiencyErrorUp(int globalBin) { return GetEfficiencyError(globalBin); };
	int FindFixBin(double x, double y = 0, double z = 0) {return effHist->FindFixBin(x,y,z); };
	TH1* GetEfficiencyHistogram() { return effHist; }
	TH1* GetPassedHistogram() { return passedHist; }
	TH1* GetTotalHistogram() { return totalHist; }
	void SetTitle(const char* title);

protected:
	TH1* passedHist;
	TH1* totalHist;
	TH1* effHist;
	Double_t totalEff;
	Double_t totalEff_error;

	ClassDef(SimpleEfficiency, 2);
};
#endif /* SIMPLEEFFICIENCY_HPP_ */
