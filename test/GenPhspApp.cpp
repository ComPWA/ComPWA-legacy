/*
 * GenPHSP.cpp
 *
 *  Created on: Oct 17, 2013
 *      Author: weidenka
 *
 *      Test app for isWithingDP function of DPKinematics. Produces uniform DP.
 */

#include <iostream>
#include <math.h>
#include <TRandom3.h>
#include <Physics/DPKinematics/DPKinematics.hpp>
#include <Physics/DPKinematics/dataPoint.hpp>
#include <TH2.h>
#include <TFile.h>
using namespace std;
int main(){
	TRandom* rand = new TRandom(54339);
	DPKinematics kin(1.869,0.0,0.497,0.493,0.493,"KS_0","K-","K+");

	TH2D* hist = new TH2D("hist","hist",100,kin.m23_sq_min,kin.m23_sq_max,100,kin.m13_sq_min,kin.m13_sq_max);

	cout<<kin.m23_sq_min<< " "<< kin.m23_sq_max<< " " << kin.m13_sq_min << " "<<kin.m13_sq_max<<endl;

	int nEvents=300000;

	int i=0;
	while( i<nEvents ){
		double m23sq = rand->Uniform(kin.m23_sq_min,kin.m23_sq_max);
		double m13sq = rand->Uniform(kin.m13_sq_min,kin.m13_sq_max);
//		double m12 = kin.getThirdVariable(sqrt(m23sq),sqrt(m13sq));
		if(!kin.isWithinDP(sqrt(m23sq),sqrt(m13sq))) continue;
		hist->Fill(m23sq,m13sq);
		i++;
	}
	TFile* tf=new TFile("test.root","recreate");
	hist->Write();
	tf->Close();
	return 1;

}




