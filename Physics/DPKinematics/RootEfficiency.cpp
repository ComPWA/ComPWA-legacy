#include "Physics/DPKinematics/RootEfficiency.hpp"
#include "Physics/DPKinematics/DataPoint.hpp"

DalitzHistEfficiency::DalitzHistEfficiency(TEfficiency* eff) : effHist(new TEfficiency(*eff)){
}
DalitzHistEfficiency::DalitzHistEfficiency(TH2D* passed, TH2D* total) : effHist(new TEfficiency(*passed, *total)){
}
DalitzHistEfficiency::DalitzHistEfficiency(const DalitzHistEfficiency&){
}
double DalitzHistEfficiency::evaluate(std::vector<double> x){
	dataPoint* point = dataPoint::instance();

	double m13sq = x[1];
	double m23sq = x[0];

	TH2D* test = (TH2D*) effHist->GetPassedHistogram();
	int globalBin = test->FindBin(m23sq,m13sq);
	return effHist->GetEfficiency(globalBin);
}
