// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

//Root header files go here
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TFile.h>

// Data Interface header files go here
#include "DIFRootReader.hpp"


using namespace std;

/************************************************************************************************/
/**
 * The main function.
 */
int DataTest(int argc, char **argv){

        DIFRootReader myReader("2Part-4vecs.root");
    unsigned int maxEvents = myReader.getNEvents();
    TLorentzVector a, b;
    double masssq;
    TH1D* bw = new TH1D("bw","inv. mass of 2 particles",1000,0.,4.);
    bw->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
    bw->GetXaxis()->CenterTitle();
    bw->GetYaxis()->SetTitle("#");
    bw->GetYaxis()->CenterTitle();

    for(unsigned int i=0; i<maxEvents; i++){
        myReader.getEvent(-1, a, b, masssq);
        bw->Fill(masssq);
    }

    TFile output("InputTest.root","RECREATE","ROOT_Tree");
    bw->Write();

    output.Write();
    output.Close();

  return 0;
}
