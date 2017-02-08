//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
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
#include "TGraph.h"
#include "TFile.h"

// Data Interface header files go here
#include "DataReader/RootReader/RootReader.hpp"
#include "DataReader/JakeReader/JakeReader.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"

//Core header files go here
#include "Core/Event.hpp"
#include "Core/Particle.hpp"

using namespace ComPWA;
using Physics::DPKinematics::DalitzKinematics;
using DataReader::RootReader::RootReader;
using DataReader::JakeReader::JakeReader;

unsigned int nBins = 400;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
	std::cout << "  ComPWA Copyright (C) 2013  Mathias Michel " << std::endl;
	std::cout << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt" << std::endl;
	std::cout << std::endl;

	std::cout << "DataIF Root 3Particles started " << std::endl << std::endl;

	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(DalitzKinematics::createInstance("jpsi", "gamma", "pi0", "pi0"));

	std::string file = "executables/Examples/DalitzFit/JPSIDATA.ACC.root";
	JakeReader myReader(file, "kin"); //return 0;
	unsigned int maxEvents = myReader.getNEvents();
	double masssq12, masssq13, masssq23;
	TH2D* bw12 = new TH2D("bw12","inv. mass-sq of particles 1&2",nBins,0.,10.,nBins,0.,10.);
	bw12->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
	bw12->GetYaxis()->SetTitle("m_{13}^{2} / GeV^{2}/c^{2}");
	bw12->GetXaxis()->CenterTitle();
	bw12->GetYaxis()->CenterTitle();
	TH2D* bw13 = new TH2D("bw13","inv. mass-sq of particles 1&3",nBins,0.,10.,nBins,0.,10.);
	bw13->GetXaxis()->SetTitle("m_{13}^{2} / GeV^{2}/c^{2}");
	bw13->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
	bw13->GetXaxis()->CenterTitle();
	bw13->GetYaxis()->CenterTitle();
	TH2D* bw23 = new TH2D("bw23","inv. mass-sq of particles 2&3",nBins,0.,10.,nBins,0.,10.);
	bw23->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}/c^{2}");
	bw23->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
	bw23->GetXaxis()->CenterTitle();
	bw23->GetYaxis()->CenterTitle();


	std::string fileB = "executables/Examples/DalitzFit/JakeJPSI_GEN.root";
	RootReader myReaderPHSP(fileB, "data");

	unsigned int maxEventsPHSP = myReaderPHSP.getNEvents();
	//double masssq12PHSP, masssq13PHSP, masssq23PHSP;
	TH2D* bw12PHSP = new TH2D("bw12GEN","inv. mass-sq of particles 1&2 GEN",nBins,0.,10.,nBins,0.,10.);
	bw12PHSP->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
	bw12PHSP->GetYaxis()->SetTitle("m_{13}^{2} / GeV^{2}/c^{2}");
	bw12PHSP->GetXaxis()->CenterTitle();
	bw12PHSP->GetYaxis()->CenterTitle();
	TH2D* bw13PHSP = new TH2D("bw13GEN","inv. mass-sq of particles 1&3 GEN",nBins,0.,10.,nBins,0.,10.);
	bw13PHSP->GetXaxis()->SetTitle("m_{13}^{2} / GeV^{2}/c^{2}");
	bw13PHSP->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
	bw13PHSP->GetXaxis()->CenterTitle();
	bw13PHSP->GetYaxis()->CenterTitle();
	TH2D* bw23PHSP = new TH2D("bw23GEN","inv. mass-sq of particles 2&3 GEN",nBins,0.,10.,nBins,0.,10.);
	bw23PHSP->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}/c^{2}");
	bw23PHSP->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
	bw23PHSP->GetXaxis()->CenterTitle();
	bw23PHSP->GetYaxis()->CenterTitle();

	double xdata[4001]; double ydata[4001];
	kin->phspContour(3,4,2000,xdata,ydata);
	TGraph* m12m13_contour = new TGraph(4001,xdata,ydata);
	m12m13_contour->SetLineColor(kRed);
	m12m13_contour->SetMarkerColor(kRed);
	m12m13_contour->SetTitle("phspContour");
	m12m13_contour->SetName("m12m13_contour");
	m12m13_contour->SetFillColor(kWhite);
	kin->phspContour(4,3,2000,xdata,ydata);
	TGraph* m13m12_contour = new TGraph(4001,xdata,ydata);
	m13m12_contour->SetLineColor(kRed);
	m13m12_contour->SetMarkerColor(kRed);
	m13m12_contour->SetTitle("phspContour");
	m13m12_contour->SetName("m13m12_contour");
	m13m12_contour->SetFillColor(kWhite);
	kin->phspContour(5,3,2000,xdata,ydata);
	TGraph* m23m12_contour = new TGraph(4001,xdata,ydata);
	m23m12_contour->SetLineColor(kRed);
	m23m12_contour->SetMarkerColor(kRed);
	m23m12_contour->SetTitle("phspContour");
	m23m12_contour->SetName("m23m13_contour");
	m23m12_contour->SetFillColor(kWhite);

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

		//m12->Fill(masssq12);
		//m23->Fill(masssq23);
		//m13->Fill(masssq13);
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

	TFile output("executables/Examples/DalitzFit/JAKEPLOTS.root","RECREATE","ROOT_Tree");
	m12m13_contour->Write("con13");
	m13m12_contour->Write("con12");
	m23m12_contour->Write("con23");
	bw12->Write("bw12");
	bw13->Write("bw13");
	bw23->Write("bw23");
	bw12PHSP->Write("phsp12");
	bw13PHSP->Write("phsp13");
	bw23PHSP->Write("phsp23");
    //m12->Write("m12");
    //m13->Write("m13");
    //m23->Write("m23");
	output.Write();
	output.Close();

	std::cout << "DataIF Root 3Particles finished " <<std::endl;

	return 0;
}
