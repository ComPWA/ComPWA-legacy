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

#include <complex>
#include <iostream>
#include "Physics/DPKinematics/SimpleEfficiency.hpp"
ClassImp(SimpleEfficiency);

SimpleEfficiency::~SimpleEfficiency(){
	delete passedHist;
	delete totalHist;
	delete effHist;
}

SimpleEfficiency::SimpleEfficiency(const SimpleEfficiency& other) : TNamed(){
//		totalEff(other.totalEff),totalEff_error(other.totalEff_error){
	((TObject&)other).Copy(*this);

	Bool_t bStatus = TH1::AddDirectoryStatus();
	TH1::AddDirectory(kFALSE);
	totalHist= (TH1*)((other.totalHist)->Clone());
	passedHist = (TH1*)((other.passedHist)->Clone());
	effHist = (TH1*)((other.effHist)->Clone());
	TH1::AddDirectory(bStatus);

	TString name = other.GetName();
	name += "_copy";
	TString title = "[copy] ";
	title += other.GetTitle();
	SetNameTitle(name,title);
}
SimpleEfficiency::SimpleEfficiency(const TH1& passed, const TH1& total){
	//we should do some consistency checks here

	Bool_t bStatus = TH1::AddDirectoryStatus();
	TH1::AddDirectory(kFALSE);
	passedHist = (TH1*) passed.Clone();
	totalHist = (TH1*) total.Clone();
	effHist = (TH1*) passed.Clone();
	effHist->Sumw2();
	effHist->Divide(totalHist);
	TH1::AddDirectory(bStatus);

	TString newName = total.GetName();
	newName += TString("_clone");
	SetName(newName);

	//	totalEff = (double)passed->GetEntries()/total->GetEntries();
	totalEff = (double)passedHist->GetSumOfWeights()/totalHist->GetSumOfWeights();
}

SimpleEfficiency::SimpleEfficiency(TString name, TString title,const TH1& passed, const TH1& total):TNamed(name,title){
	//we should do some consistency checks here

	Bool_t bStatus = TH1::AddDirectoryStatus();
	TH1::AddDirectory(kFALSE);
	passedHist = (TH1*) passed.Clone(title+" (passed)");
	totalHist = (TH1*) total.Clone(title+" (total)");
	effHist = (TH1*) passed.Clone(title+" (eff)");
	effHist->Sumw2();
	effHist->Divide(totalHist);
	TH1::AddDirectory(bStatus);

	//	totalEff = (double)passed->GetEntries()/total->GetEntries();
	totalEff = (double)passedHist->GetSumOfWeights()/totalHist->GetSumOfWeights();
	totalEff_error = totalEff*sqrt(1/passedHist->GetSumOfWeights()+1/totalHist->GetSumOfWeights());
};
double SimpleEfficiency::GetEfficiency(int globalBin){
	return (double) effHist->GetBinContent(globalBin);
}
double SimpleEfficiency::GetEfficiencyError(int globalBin){
	unsigned int N = totalHist->GetSumOfWeights();
	double eff = GetEfficiency(globalBin);
	double err = 1;
	return err;
}
void SimpleEfficiency::SetTitle(const char* title){
	//setting the titles (looking for the first semicolon and insert the tokens there)
	TString title_passed = title;
	TString title_total = title;
	Ssiz_t pos = title_passed.First(";");
	if (pos != kNPOS) {
		title_passed.Insert(pos," (passed)");
		title_total.Insert(pos," (total)");
	}
	else {
		title_passed.Append(" (passed)");
		title_total.Append(" (total)");
	}
	passedHist->SetTitle(title_passed);
	totalHist->SetTitle(title_total);
	effHist->SetTitle(TString(title)+";efficiency");

	// strip (total) for the TEfficiency title
	// HIstogram SetTitle has already stripped the axis
	TString teffTitle = totalHist->GetTitle();
	teffTitle.ReplaceAll(" (total)","");
	TNamed::SetTitle(teffTitle);
	return;
}
