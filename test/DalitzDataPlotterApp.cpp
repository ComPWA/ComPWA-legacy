//! Plotter for three particle final states.
/*!
 * @file DalitzDataPlotterApp.cpp
 * This tiny application uses the RootReader to read a file with three final
 * particles and plots the invariant masses of the two particle subsystems.
 * The resulting Dalitz-plots are written as histograms to another root file.
*/

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

//Root header files go here
#include "TLorentzVector.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"

// Data Interface header files go here
#include "DataReader/RootReader/RootReader.hpp"

//Core header files go here
#include "Core/Event.hpp"
#include "Core/Particle.hpp"

using namespace std;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){

	cout << "DataIF Root 3Particles started " << endl << endl;
	double max_range=2.;
	double min_range=.9;
	int nBins=300;
	string file = "test/3Part-4vecs.root";
	RootReader myReader(file, false,"data");
	unsigned int maxEvents = myReader.getNEvents();
	double masssq12, masssq13, masssq23;
	TH2D* bw1213 = new TH2D("bw1213","m12 versus m13",nBins,min_range,max_range,nBins,min_range,max_range);
	bw1213->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
	bw1213->GetXaxis()->CenterTitle();
	bw1213->GetYaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
	bw1213->GetYaxis()->CenterTitle();
	TH2D* bw1312 = new TH2D("bw1312","m13 versus m12",nBins,min_range,max_range,nBins,min_range,max_range);
	bw1312->GetXaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
	bw1312->GetXaxis()->CenterTitle();
	bw1312->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
	bw1312->GetYaxis()->CenterTitle();
	TH2D* bw2313 = new TH2D("bw2313","m23 versus m13",nBins,min_range,max_range,nBins,min_range,max_range);
	bw2313->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}");
	bw2313->GetXaxis()->CenterTitle();
	bw2313->GetYaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
	bw2313->GetYaxis()->CenterTitle();
	TH2D* bw2312 = new TH2D("bw2312","m23 versus m12",nBins,min_range,max_range,nBins,min_range,max_range);
	bw2312->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}");
	bw2312->GetXaxis()->CenterTitle();
	bw2312->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
	bw2312->GetYaxis()->CenterTitle();

    RootReader myReaderPHSP(file, false,"mc");
    unsigned int maxEventsPHSP = myReaderPHSP.getNEvents();
    //double masssq12PHSP, masssq13PHSP, masssq23PHSP;
    TH2D* bw12PHSP = new TH2D("bw12PHSP","inv. mass-sq of particles 1&2 PHSP",nBins,min_range,max_range,nBins,min_range,max_range);
    bw12PHSP->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
    bw12PHSP->GetXaxis()->CenterTitle();
    bw12PHSP->GetYaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
    bw12PHSP->GetYaxis()->CenterTitle();
    TH2D* bw13PHSP = new TH2D("bw13PHSP","inv. mass-sq of particles 1&3 PHSP",nBins,min_range,max_range,nBins,min_range,max_range);
    bw13PHSP->GetXaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
    bw13PHSP->GetXaxis()->CenterTitle();
    bw13PHSP->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
    bw13PHSP->GetYaxis()->CenterTitle();
    TH2D* bw23PHSP = new TH2D("bw23PHSP","inv. mass-sq of particles 2&3 PHSP",nBins,min_range,max_range,nBins,min_range,max_range);
    bw23PHSP->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}");
    bw23PHSP->GetXaxis()->CenterTitle();
    bw23PHSP->GetYaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
    bw23PHSP->GetYaxis()->CenterTitle();

    for(unsigned int i = 0; i < maxEvents; i++){
        Event event(myReader.getEvent(i));

    	//myReader.getEvent(-1, a, b, masssq);
    	//if(!myReader.getEvent(i, event)) continue; TODO: try exception
    	if(!event.getNParticles() == 3) continue;
    	//if(!event) continue;
    	//cout << "Event: \t" << i << "\t NParticles: \t" << event.getNParticles() << endl;
    	const Particle &a(event.getParticle(0));
    	const Particle &b(event.getParticle(1));
    	const Particle &c(event.getParticle(2));
    	masssq12 = pow(a.E+b.E,2) - pow(a.px+b.px ,2) - pow(a.py+b.py ,2) - pow(a.pz+b.pz ,2);
        masssq13 = pow(a.E+c.E,2) - pow(a.px+c.px ,2) - pow(a.py+c.py ,2) - pow(a.pz+c.pz ,2);
        masssq23 = pow(b.E+c.E,2) - pow(b.px+c.px ,2) - pow(b.py+c.py ,2) - pow(b.pz+c.pz ,2);

    	bw1213->Fill(masssq12,masssq13);
        bw1312->Fill(masssq13,masssq12);
        bw2313->Fill(masssq23,masssq13);
        bw2312->Fill(masssq23,masssq12);
    }

    for(unsigned int i = 0; i < maxEventsPHSP; i++){
        Event event(myReaderPHSP.getEvent(i));

        //myReader.getEvent(-1, a, b, masssq);
        //if(!myReader.getEvent(i, event)) continue; TODO: try exception
        if(!event.getNParticles() == 3) continue;
        //if(!event) continue;
        //cout << "Event: \t" << i << "\t NParticles: \t" << event.getNParticles() << endl;
        const Particle &a(event.getParticle(0));
        const Particle &b(event.getParticle(1));
        const Particle &c(event.getParticle(2));
        masssq12 = pow(a.E+b.E,2) - pow(a.px+b.px ,2) - pow(a.py+b.py ,2) - pow(a.pz+b.pz ,2);
        masssq13 = pow(a.E+c.E,2) - pow(a.px+c.px ,2) - pow(a.py+c.py ,2) - pow(a.pz+c.pz ,2);
        masssq23 = pow(b.E+c.E,2) - pow(b.px+c.px ,2) - pow(b.py+c.py ,2) - pow(b.pz+c.pz ,2);

        bw12PHSP->Fill(masssq12,masssq13);
        bw13PHSP->Fill(masssq13,masssq12);
        bw23PHSP->Fill(masssq23,masssq13);
    }

    TFile output("out.root","RECREATE","ROOT_Tree");
    bw1213->Write();
    bw1312->Write();
    bw2313->Write();
    bw2312->Write();
    bw12PHSP->Write();
    bw13PHSP->Write();
    bw23PHSP->Write();
    TCanvas* c1 = new TCanvas("c1","dalitz distributions",200,10,1400,700);
    c1->Divide(3,2);
    c1->cd(1); bw2312->Draw("COLZ");
    c1->cd(2); bw2313->Draw("COLZ");
    c1->cd(3); bw1213->Draw("COLZ");
    c1->cd(4); bw2312->ProjectionX()->DrawCopy();
    c1->cd(5); bw2312->ProjectionY()->DrawCopy();
    c1->cd(6); bw1312->ProjectionX()->DrawCopy();
    c1->Write();
    output.Write();
    output.Close();

    cout << "DataIF Root 3Particles finished " <<endl;

  return 0;
}
