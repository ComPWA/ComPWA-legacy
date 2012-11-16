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
#include "DIFRootReader.hpp"

//Core header files go here
#include "PWAEvent.hpp"
#include "PWAParticle.hpp"

using namespace std;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){

    cout << "DataIF Root 2Particles started " << endl << endl;

    DIFRootReader myReader("test/2Part-4vecs.root");
    unsigned int maxEvents = myReader.getNEvents();
    double masssq;
    TH1D* bw = new TH1D("bw","inv. mass of 2 particles",1000,0.,4.);
    bw->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
    bw->GetXaxis()->CenterTitle();
    bw->GetYaxis()->SetTitle("#");
    bw->GetYaxis()->CenterTitle();

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

    	bw->Fill(masssq);
    }

    TFile output("test/InputTest.root","RECREATE","ROOT_Tree");
    bw->Write();

    output.Write();
    output.Close();

    cout << "DataIF Root 2Particles finished " <<endl;

  return 0;
}