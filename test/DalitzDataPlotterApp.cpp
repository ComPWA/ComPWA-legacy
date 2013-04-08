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
#include "Core/PWAEvent.hpp"
#include "Core/PWAParticle.hpp"

using namespace std;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){

    cout << "DataIF Root 3Particles started " << endl << endl;

    string file = "test/3Part-4vecs.root";
    RootReader myReader(file, false);
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

    for(unsigned int i = 0; i < maxEvents; i++){
        PWAEvent event(myReader.getEvent(i));
        PWAParticle a, b, c;

    	//myReader.getEvent(-1, a, b, masssq);
    	//if(!myReader.getEvent(i, event)) continue; TODO: try exception
    	if(!event.getNParticles() == 3) continue;
    	//if(!event) continue;
    	//cout << "Event: \t" << i << "\t NParticles: \t" << event.getNParticles() << endl;
    	event.getParticle(0,a);
    	event.getParticle(1,b);
        event.getParticle(2,c);
    	masssq12 = pow(a.getE()+b.getE(),2) - pow(a.getPx()+b.getPx() ,2) - pow(a.getPy()+b.getPy() ,2) - pow(a.getPz()+b.getPz() ,2);
        masssq13 = pow(a.getE()+c.getE(),2) - pow(a.getPx()+c.getPx() ,2) - pow(a.getPy()+c.getPy() ,2) - pow(a.getPz()+c.getPz() ,2);
        masssq23 = pow(b.getE()+c.getE(),2) - pow(b.getPx()+c.getPx() ,2) - pow(b.getPy()+c.getPy() ,2) - pow(b.getPz()+c.getPz() ,2);

    	bw12->Fill(masssq12,masssq13);
        bw13->Fill(masssq13,masssq12);
        bw23->Fill(masssq23,masssq12);
    }

    TFile output("test/DalitzJPSI.root","RECREATE","ROOT_Tree");
    bw12->Write();
    bw13->Write();
    bw23->Write();
    output.Write();
    output.Close();

    cout << "DataIF Root 3Particles finished " <<endl;

  return 0;
}
