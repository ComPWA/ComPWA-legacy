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
#include "Physics/DPKinematics/RootEfficiency.hpp"
#include "Core/DataPoint.hpp"

DalitzHistEfficiency::DalitzHistEfficiency(TEfficiency* eff) : effHist(new TEfficiency(*eff)){
	BOOST_LOG_TRIVIAL(debug) << "DalitzHistEfficiency: creating efficiency from existing TEfficiency object!";
}
DalitzHistEfficiency::DalitzHistEfficiency(TH2D* passed, TH2D* total) : effHist(new TEfficiency(*passed, *total)){
	BOOST_LOG_TRIVIAL(debug) << "DalitzHistEfficiency: creating efficiency from two TH2D objects!";
}
DalitzHistEfficiency::DalitzHistEfficiency(const DalitzHistEfficiency&){
}
double DalitzHistEfficiency::evaluate(std::vector<double> x){
	dataPoint point; point.setVal("m23sq",x[0]); point.setVal("m13sq",x[1]);
//	double m13sq = x[1];
//	double m23sq = x[0];

//	TH2D* test = (TH2D*) effHist->GetPassedHistogram();
//	int globalBin = test->FindBin(m23sq,m13sq);
//	return effHist->GetEfficiency(globalBin);
	return evaluate(point);
}
double DalitzHistEfficiency::evaluate(dataPoint& point){
	double m13sq = point.getVal("m13sq");
	double m23sq = point.getVal("m23sq");

	TH2D* test = (TH2D*) effHist->GetPassedHistogram();
	int globalBin = test->FindBin(m23sq,m13sq);
	return effHist->GetEfficiency(globalBin);
}

double DalitzAngleHistEfficiency::evaluate(dataPoint& point){
	double m23sq = point.getVal("m23sq");
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double angle = kin->calcHelicityAngle(5,point);

	TH2D* test = (TH2D*) effHist->GetPassedHistogram();
	int globalBin = test->FindBin(m23sq,angle);
	return effHist->GetEfficiency(globalBin);
}
