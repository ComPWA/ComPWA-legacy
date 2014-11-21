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


#include "Core/DataPoint.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"
#include "Physics/DPKinematics/SimpleAngleEfficiency.hpp"


//SimpleAngleEfficiency::SimpleAngleEfficiency(SimpleEfficiency* eff) : effHist(new SimpleEfficiency(*eff)){
SimpleAngleEfficiency::SimpleAngleEfficiency(SimpleEfficiency* eff) {
	effHist = std::shared_ptr<SimpleEfficiency>(new SimpleEfficiency(*eff));
	BOOST_LOG_TRIVIAL(debug) << "SimpleAngleEfficiency: creating efficiency from existing SimpleEfficiency object!";
}
SimpleAngleEfficiency::SimpleAngleEfficiency(TH1* passed, TH1* total) : effHist(new SimpleEfficiency("angleEff","angleEff",*passed, *total)){
	BOOST_LOG_TRIVIAL(debug) << "SimpleAngleEfficiency: creating efficiency from two histograms objects!";
}
SimpleAngleEfficiency::SimpleAngleEfficiency(const SimpleAngleEfficiency& other){
	effHist = other.effHist;
}
double SimpleAngleEfficiency::evaluate(std::vector<double> x){
	//assume that x[0]=m23sq and x[1]=m13sq
	dataPoint point; point.setVal(0,x[0]); point.setVal(1,x[1]);
	return evaluate(point);
}
double SimpleAngleEfficiency::evaluate(dataPoint& point){
	double m23sq = point.getVal(0);
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double angle = kin->helicityAngle(5,point);

	TH1* test = (TH1*) effHist->GetPassedHistogram();
	int globalBin = test->FindBin(m23sq,angle);
	double eff = effHist->GetEfficiency(globalBin);
	return eff;
}
