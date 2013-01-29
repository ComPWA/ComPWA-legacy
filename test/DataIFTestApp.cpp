//! Test-Application of the Root Data-IF.
/*!
 * @file DataIFTestApp.cpp
 * This tiny application tests the Root-format Data-IF. It reads a root file,
 * plots the invariant mass of the first two particles and writes the result as
 * well to another root file.
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
#include "TH1D.h"
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

    cout << "DataIF Root 2Particles started " << endl << endl;

    string file = "test/2Part-4vecs.root";
    RootReader myReader(file, false);
    unsigned int maxEvents = myReader.getNEvents();
    double masssq;
    TH1D* bw = new TH1D("bw","inv. mass of 2 particles",1000,0.,2.4);
    bw->GetXaxis()->SetTitle("m_{12} / GeV");
    bw->GetXaxis()->CenterTitle();
    bw->GetYaxis()->SetTitle("#");
    bw->GetYaxis()->CenterTitle();
    TH1D* bw2 = new TH1D("bw2","inv. mass-sq of 2 particles",1000,0.,4.);
    bw2->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
    bw2->GetXaxis()->CenterTitle();
    bw2->GetYaxis()->SetTitle("#");
    bw2->GetYaxis()->CenterTitle();

    for(unsigned int i = 0; i < maxEvents; i++){
        PWAEvent event;
        PWAParticle a, b;

    	//myReader.getEvent(-1, a, b, masssq);
    	if(!myReader.getEvent(i, event)) continue;
    	//if(!event) continue;
    	//cout << "Event: \t" << i << "\t NParticles: \t" << event.getNParticles() << endl;
    	event.getParticle(0,a);
    	event.getParticle(1,b);
    	masssq = pow(a.getE()+b.getE(),2) - pow(a.getPx()+b.getPx() ,2) - pow(a.getPy()+b.getPy() ,2) - pow(a.getPz()+b.getPz() ,2);

        bw->Fill(sqrt(masssq));
    	bw2->Fill(masssq);
    }

    TFile output("test/InputTest.root","RECREATE","ROOT_Tree");
    bw->Write();
    bw2->Write();
    output.Write();
    output.Close();

    cout << "DataIF Root 2Particles finished " <<endl;

  return 0;
}
