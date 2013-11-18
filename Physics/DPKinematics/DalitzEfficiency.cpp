#include "Physics/DPKinematics/DalitzEfficiency.hpp"
#include "Physics/DPKinematics/DataPoint.hpp"

DalitzEfficiency::DalitzEfficiency(){
}
DalitzEfficiency::DalitzEfficiency(const DalitzEfficiency&, const char*){
}
DalitzEfficiency::~DalitzEfficiency(){
}

DalitzHistEfficiency::DalitzHistEfficiency(TEfficiency* eff) : effHist(new TEfficiency(*eff)){
}
DalitzHistEfficiency::DalitzHistEfficiency(TH2D* passed, TH2D* total) : effHist(new TEfficiency(*passed, *total)){
}
DalitzHistEfficiency::DalitzHistEfficiency(const DalitzEfficiency&){
}
double DalitzHistEfficiency::evaluate(){
	dataPoint* point = dataPoint::instance();

	double m13sq = point->getMsq(4);
	double m23sq = point->getMsq(5);

	TH2D* test = (TH2D*) effHist->GetPassedHistogram();
	int globalBin = test->FindBin(m23sq,m13sq);
	return effHist->GetEfficiency(globalBin);
}
