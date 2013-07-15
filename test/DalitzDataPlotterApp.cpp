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

    string file = "test/3Part-4vecs.root";
    RootReader myReader(file, false,"data");
    unsigned int maxEvents = myReader.getNEvents();
    double masssq12, masssq13, masssq23;
    TH2D* bw12 = new TH2D("bw12","inv. mass-sq of particles 1&2",1000,0.,10.,1000,0.,10.);
    bw12->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
    bw12->GetXaxis()->CenterTitle();
    bw12->GetYaxis()->SetTitle("#");
    bw12->GetYaxis()->CenterTitle();
    TH2D* bw13 = new TH2D("bw13","inv. mass-sq of particles 1&3",1000,0.,10.,1000,0.,10.);
    bw13->GetXaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
    bw13->GetXaxis()->CenterTitle();
    bw13->GetYaxis()->SetTitle("#");
    bw13->GetYaxis()->CenterTitle();
    TH2D* bw23 = new TH2D("bw23","inv. mass-sq of particles 2&3",1000,0.,10.,1000,0.,10.);
    bw23->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}");
    bw23->GetXaxis()->CenterTitle();
    bw23->GetYaxis()->SetTitle("#");
    bw23->GetYaxis()->CenterTitle();

    RootReader myReaderPHSP(file, false,"mc");
    unsigned int maxEventsPHSP = myReaderPHSP.getNEvents();
    //double masssq12PHSP, masssq13PHSP, masssq23PHSP;
    TH2D* bw12PHSP = new TH2D("bw12PHSP","inv. mass-sq of particles 1&2 PHSP",1000,0.,10.,1000,0.,10.);
    bw12PHSP->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
    bw12PHSP->GetXaxis()->CenterTitle();
    bw12PHSP->GetYaxis()->SetTitle("#");
    bw12PHSP->GetYaxis()->CenterTitle();
    TH2D* bw13PHSP = new TH2D("bw13PHSP","inv. mass-sq of particles 1&3 PHSP",1000,0.,10.,1000,0.,10.);
    bw13PHSP->GetXaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
    bw13PHSP->GetXaxis()->CenterTitle();
    bw13PHSP->GetYaxis()->SetTitle("#");
    bw13PHSP->GetYaxis()->CenterTitle();
    TH2D* bw23PHSP = new TH2D("bw23PHSP","inv. mass-sq of particles 2&3 PHSP",1000,0.,10.,1000,0.,10.);
    bw23PHSP->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}");
    bw23PHSP->GetXaxis()->CenterTitle();
    bw23PHSP->GetYaxis()->SetTitle("#");
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

    	bw12->Fill(masssq12,masssq13);
        bw13->Fill(masssq13,masssq12);
        bw23->Fill(masssq23,masssq12);
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
        bw23PHSP->Fill(masssq23,masssq12);
    }

    TFile output("test/DalitzJPSI.root","RECREATE","ROOT_Tree");
    bw12->Write();
    bw13->Write();
    bw23->Write();
    bw12PHSP->Write();
    bw13PHSP->Write();
    bw23PHSP->Write();
    output.Write();
    output.Close();

    cout << "DataIF Root 3Particles finished " <<endl;

  return 0;
}
