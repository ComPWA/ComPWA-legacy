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
//! Application to generate J/Psi -> g pi pi.
/*!
 * @file GenDalitzApp.cpp
 * This application uses the Breit-Wigner-Amplitude-Sum module and a
 * phase-space generator to generate a file with J/Psi -> g pi pi events.
 */

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

// Boost header files go here
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
using namespace boost::log;

// Root header files go here
#include "TMath.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TGenPhaseSpace.h"
#include "TFile.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TRandom3.h"

// Physics Interface header files go here
#include "Core/PhysConst.hpp"
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Efficiency.hpp"

using namespace std;

const unsigned int MaxEvents = 10;

//constants

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
	boost::log::core::get()->set_filter(trivial::severity >= trivial::info); //setting log level
	BOOST_LOG_TRIVIAL(info)<< "  ComPWA Copyright (C) 2013  Mathias Michel ";
	BOOST_LOG_TRIVIAL(info)<< "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt";
	BOOST_LOG_TRIVIAL(info)<<"";

	unsigned int i=0, mc=0;
	TRandom3 rando;

	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(DalitzKinematics::createInstance("J/psi","gamma","pi0","pi0"));
//	static dataPoint* point = dataPoint::instance();


	//DPKinematics kin("J/psi","gamma","pi0","pi0");
	//DPKinematics kin("D0","gamma","K-","K+");
	//static dataPoint* point = dataPoint::instance(kin);

	/*const double M = kin.getMass("D0"); // GeV/c² (J/psi+)
	const double Br = 0.000093; // GeV/c² (width)
	const double m1 = kin.getMass("gamma"); // GeV/c² (gamma)
	const double m2 = kin.getMass("K-"); // GeV/c² (pi)
	const double m3 = kin.getMass("K+"); // GeV/c² (pi)
	//const double c = 299792458.; // m/s
	const double PI = PhysConst::instance()->getConstValue("Pi");*/

//	const Double_t M = 3.096916; // GeV/c² (J/psi+)
//	const Double_t Br = 0.000093; // GeV/c² (width)
//	const Double_t m1 = 0.; // GeV/c² (gamma)
//	const Double_t m2 = 0.139570; // GeV/c² (pi)
//	const Double_t m3 = 0.139570; // GeV/c² (pi)
	//const Double_t c = 299792458.; // m/s
	/*const Double_t PI = 3.14159; // m/s

  //load resonances
  //DPKinematics kin(M,Br,m1,m2,m3,"gamma","pi0","pi0");
  //dataPoint::instance(kin);
  std::string resoFile="test/JPSI_ypipi.xml";
  AmplitudeSetup ini(resoFile);
  cout << "loaded file " << ini.getFileName() << " with " << ini.getResonances().size() << " resonances:" << endl;
  for(std::vector<Resonance>::iterator reso=ini.getResonances().begin(); reso!=ini.getResonances().end(); reso++){
    cout << endl << "Resonance " << (*reso).m_name << endl;
    cout << "Mass =  " << (*reso).m_mass << " with range " << (*reso).m_mass_min << " to " << (*reso).m_mass_max << endl;
    cout << "Width = " << (*reso).m_width << " with range " << (*reso).m_width_min << " to " << (*reso).m_width_max << endl;
    cout << "Spin =  " << (*reso).m_spin << " m = " << (*reso).m_m << " n = " << (*reso).m_n << endl;
    cout << "Strength =  " << (*reso).m_strength << " Phase = " << (*reso).m_phase << endl;
    cout << "Breakupmomentum =  " << (*reso).m_mesonRadius<< endl;
    cout << "DaughterA =  " << (*reso).m_daugtherA << " DaughterB = " << (*reso).m_daugtherB << endl;
  }
  cout << endl << endl;

  //Simple Breit-Wigner Physics-Module setup
  AmpSumIntensity testBW(kin, ini, MaxEvents, AmpSumIntensity::normalizationStyle::none);
  cout << testBW.printAmps() << endl;

 // return 0;

  ParameterList minPar;
  testBW.fillStartParVec(minPar);
  cout << minPar << endl;
 // minPar.AddParameter(DoubleParameter(1.5,0.5,2.5,0.1));

  //Output File setup
  TFile output("test/3Part-4vecs.root","recreate");
  output.SetCompressionLevel(1); //try level 2 also

  TTree fTree ("data","Dalitz-Gen");
  TClonesArray *fEvt = new TClonesArray("TParticle");
  //TClonesArray &ar = *fEvt;
  fTree.Branch("Particles",&fEvt);

  TTree fTreePHSP ("mc","Dalitz-Gen-PHSP");
  TClonesArray *fEvtPHSP = new TClonesArray("TParticle");
  //TClonesArray &ar = *fEvt;
  fTreePHSP.Branch("Particles",&fEvtPHSP);

  //Generation
  TLorentzVector W(0.0, 0.0, 0.0, M);//= beam + target;

  //(Momentum, Energy units are Gev/C, GeV)
  Double_t masses[3] = { m1, m2, m2} ;

  TGenPhaseSpace event;
  event.SetDecay(W, 3, masses);

  TLorentzVector *pGamma,*pPip,*pPim,pPm23,pPm13,pPm12;
  double weight, m23sq, m13sq, m12sq, maxTest=0;
  cout << "Einschwingen" << endl;
  for(unsigned int schwing=0; schwing<10*MaxEvents; schwing++){
      weight = event.Generate();

      pGamma = event.GetDecay(0);
      pPip    = event.GetDecay(1);
      pPim    = event.GetDecay(2);

	  pPm23 = *pPim + *pPip;
	  pPm13 = *pGamma + *pPim;
	  pPm12 = *pGamma + *pPip;

      m23sq=pPm23.M2(); m13sq=pPm13.M2(); m12sq=pPm12.M2();

      //		m12sq = kin.getThirdVariableSq(m23sq,m13sq);
      		point->setMsq(3,m12sq); point->setMsq(4,m13sq); point->setMsq(5,m23sq);
      //		m12sq=M*M+m1*m1+m2*m2+m3*m3-m13sq-m23sq;
      		if( abs(m12sq-kin.getThirdVariableSq(m23sq,m13sq))>0.01 ){
      			std::cout<<m12sq<<" "<<kin.getThirdVariableSq(m23sq,m13sq)<<std::endl;
      			std::cout<<"   " <<m23sq<<" "<<m13sq<<" "<<m12sq<<std::endl;
      		}

      //call physics module
      vector<double> x;
      x.push_back(sqrt(m23sq));
      x.push_back(sqrt(m13sq));
      x.push_back(sqrt(m12sq));
      ParameterList intensL = testBW.intensity(x, minPar);
      double AMPpdf = intensL.GetDoubleParameter(0)->GetValue();
      //double AMPpdf = testBW.intensity(x, minPar);


      //mb.setVal(m13);
      //double m13pdf = totAmp13.getVal();//fun_combi2->Eval(m13);
      if(maxTest<(weight*AMPpdf))
        maxTest=(weight*AMPpdf);

  }

  maxTest*=1.1;
  int outCnt=0, maxCnt=MaxEvents/100;
    cout << "Start generation of y pi0 pi0 Dalitz" << endl;
    do{
        weight = event.Generate();

        pGamma = event.GetDecay(0);
        pPip    = event.GetDecay(1);
        pPim    = event.GetDecay(2);

        pPm23 = *pPim + *pPip;
        pPm13 = *pGamma + *pPim;
        pPm12 = *pGamma + *pPip;

        m23sq=pPm23.M2(); m13sq=pPm13.M2(); m12sq=pPm12.M2();

        //		m12sq = kin.getThirdVariableSq(m23sq,m13sq);
        		point->setMsq(3,m12sq); point->setMsq(4,m13sq); point->setMsq(5,m23sq);
        //		m12sq=M*M+m1*m1+m2*m2+m3*m3-m13sq-m23sq;
          		if( abs(m12sq-kin.getThirdVariableSq(m23sq,m13sq))>0.01 ){
          			std::cout<<m12sq<<" "<<kin.getThirdVariableSq(m23sq,m13sq)<<std::endl;
          			std::cout<<"   " <<m23sq<<" "<<m13sq<<" "<<m12sq<<std::endl;
          		}
        TParticle fparticleGam(22,1,0,0,0,0,*pGamma,W);
        TParticle fparticlePip(211,1,0,0,0,0,*pPip,W);
        TParticle fparticlePim(-211,1,0,0,0,0,*pPim,W);

        //call physics module
        vector<double> x;
        x.push_back(sqrt(m23sq));
        x.push_back(sqrt(m13sq));
        x.push_back(sqrt(m12sq));
        ParameterList intensL = testBW.intensity(x, minPar);
        double AMPpdf = intensL.GetDoubleParameter(0)->GetValue();
        //double AMPpdf = testBW.intensity(x, minPar);

        double test = rando.Uniform(0,maxTest);
        double testmc = rando.Uniform(0,1.);

        //mb.setVal(m13);
        //double m13pdf = totAmp13.getVal();//fun_combi2->Eval(m13);
        if(maxTest<(weight*AMPpdf))
          cout << "Einschwingen zu kurz!" << endl;
        if(i<MaxEvents && test<(weight*AMPpdf)){
          if(outCnt==maxCnt){
        	outCnt=0;
            cout << (i/(double)MaxEvents*100.) << "% : " << test << " " << (weight*AMPpdf) << endl;
          }
          outCnt++;
          i++;
          new((*fEvt)[0]) TParticle(fparticleGam);
          new((*fEvt)[1]) TParticle(fparticlePip);
          new((*fEvt)[2]) TParticle(fparticlePim);

          fTree.Fill();
        }

        if(mc<MaxEvents && testmc<weight){
          mc++;
          new((*fEvtPHSP)[0]) TParticle(fparticleGam);
          new((*fEvtPHSP)[1]) TParticle(fparticlePip);
          new((*fEvtPHSP)[2]) TParticle(fparticlePim);

          fTreePHSP.Fill();
        }
    }while(i<MaxEvents || mc<MaxEvents);
    cout << "100%! Write Data" << endl;

  fTree.Print();
  fTree.Write();
  fTreePHSP.Write();
  output.Close();

  cout << testBW.printAmps() << endl;
  cout << "Done ... " << maxTest << endl << endl;

  return 0;*/

//	const Double_t PI = 3.14159; // m/s

	//load resonances
	//DPKinematics kin(M,Br,m1,m2,m3,"gamma","pi0","pi0");
	//dataPoint::instance(kin);

	const char* pPath = getenv("COMPWA_DIR");
	std::string path = "";
	try{
	  path = std::string(pPath);
	}catch(std::logic_error){
	  BOOST_LOG_TRIVIAL(error)<<"Environment Variable COMPWA_DIR not set?"<<std::endl;
	}
	std::string resoFile=path+"/test/JPSI_ypipi.xml";
	AmplitudeSetup ini(resoFile);
	cout << "loaded file " << ini.getFileName() << " with " << ini.getResonances().size() << " resonances:" << endl;
	for(std::vector<Resonance>::iterator reso=ini.getResonances().begin(); reso!=ini.getResonances().end(); reso++){
		cout << endl << "Resonance " << (*reso).m_name << endl;
		cout << "Mass =  " << (*reso).m_mass << " with range " << (*reso).m_mass_min << " to " << (*reso).m_mass_max << endl;
		cout << "Width = " << (*reso).m_width << " with range " << (*reso).m_width_min << " to " << (*reso).m_width_max << endl;
		cout << "Spin =  " << (*reso).m_spin << " m = " << (*reso).m_m << " n = " << (*reso).m_n << endl;
		cout << "Strength =  " << (*reso).m_strength << " Phase = " << (*reso).m_phase << endl;
		cout << "mesonRadius=  " << (*reso).m_mesonRadius<< endl;
		cout << "DaughterA =  " << (*reso).m_daugtherA << " DaughterB = " << (*reso).m_daugtherB << endl;
	}
	cout << endl << endl;

	//Simple Breit-Wigner Physics-Module setup
	AmpSumIntensity testBW(ini, AmpSumIntensity::normStyle::none, std::shared_ptr<Efficiency>(new UnitEfficiency()), MaxEvents);
	 // std::shared_ptr<Amplitude> amps(new AmpSumIntensity(ini, AmpSumIntensity::normStyle::one, std::shared_ptr<Efficiency>(new UnitEfficiency()), myReader->getNEvents()));

	testBW.printAmps();

	ParameterList minPar;
	testBW.fillStartParVec(minPar);
	cout << minPar << endl;
	// minPar.AddParameter(DoubleParameter(1.5,0.5,2.5,0.1));

	//Output File setup
	TFile output("test/TEST.root","recreate");
	//TFile output("test/3Part-4vecs.root","recreate");
	output.SetCompressionLevel(1); //try level 2 also

	TTree fTree ("data","Dalitz-Gen");
	TClonesArray *fEvt = new TClonesArray("TParticle");
	//TClonesArray &ar = *fEvt;
	fTree.Branch("Particles",&fEvt);

	TTree fTreePHSP ("mc","Dalitz-Gen-PHSP");
	TClonesArray *fEvtPHSP = new TClonesArray("TParticle");
	//TClonesArray &ar = *fEvt;
	fTreePHSP.Branch("Particles",&fEvtPHSP);

	//Generation
//	TLorentzVector W(0.0, 0.0, 0.0, M);//= beam + target;
	TLorentzVector W(0.0, 0.0, 0.0, kin->M);//= beam + target;

	//(Momentum, Energy units are Gev/C, GeV)
	Double_t masses[3] = { kin->m1, kin->m2, kin->m2} ;

	TGenPhaseSpace event;
	event.SetDecay(W, 3, masses);

	TLorentzVector *pGamma,*pPip,*pPim,pPm23,pPm13,pPm12;
	double weight, m23sq, m13sq, m12sq, maxTest=0, m12sqtest;
	cout << "Einschwingen" << endl;
	for(unsigned int schwing=0; schwing<10*MaxEvents; schwing++){
		weight = event.Generate();

		pGamma = event.GetDecay(0);
		pPip    = event.GetDecay(1);
		pPim    = event.GetDecay(2);

		pPm23 = *pPim + *pPip;
		pPm13 = *pGamma + *pPim;
		pPm12 = *pGamma + *pPip;

		m23sq=pPm23.M2(); m13sq=pPm13.M2(); m12sq=pPm12.M2();

	dataPoint dataP; dataP.setVal("m23sq",m23sq);	dataP.setVal("m13sq",m13sq);
//		point->setMsq(3,m12sq); point->setMsq(4,m13sq); point->setMsq(5,m23sq);
		//		m12sq=M*M+m1*m1+m2*m2+m3*m3-m13sq-m23sq;
	    m12sqtest = kin->getThirdVariableSq(m23sq,m13sq);
		if( (m12sq-m12sqtest)>0.01 || (m12sqtest-m12sq)>0.01 ){
			std::cout<<m12sq<<" "<<m12sqtest<<std::endl;
			std::cout<<"   " <<m23sq<<" "<<m13sq<<" "<<m12sq<<std::endl;
		}

		//call physics module
//		vector<double> x;
//		x.push_back(sqrt(m23sq));
//		x.push_back(sqrt(m13sq));
//		x.push_back(sqrt(m12sq));
		//ParameterList intensL = testBW.intensity(dataP);
		double AMPpdf = testBW.intensity(dataP).GetParameterValue(0);
		//double AMPpdf = testBW.intensity(x, minPar);


		//mb.setVal(m13);
		//double m13pdf = totAmp13.getVal();//fun_combi2->Eval(m13);
		if(maxTest<(weight*AMPpdf))
			maxTest=(weight*AMPpdf);

	}

	maxTest*=1.1;
	int outCnt=0, maxCnt=MaxEvents/100;
	cout << "Start generation of y pi0 pi0 Dalitz" << endl;
	do{
		weight = event.Generate();

		pGamma = event.GetDecay(0);
		pPip    = event.GetDecay(1);
		pPim    = event.GetDecay(2);

		pPm23 = *pPim + *pPip;
		pPm13 = *pGamma + *pPim;
		pPm12 = *pGamma + *pPip;

		m23sq=pPm23.M2(); m13sq=pPm13.M2(); m12sq=pPm12.M2();

	dataPoint dataP; dataP.setVal("m23sq",m23sq);	dataP.setVal("m13sq",m13sq);
		//		m12sq = kin.getThirdVariableSq(m23sq,m13sq);
//		point->setMsq(3,m12sq); point->setMsq(4,m13sq); point->setMsq(5,m23sq);
		//		m12sq=M*M+m1*m1+m2*m2+m3*m3-m13sq-m23sq;
		if( abs(m12sq-kin->getThirdVariableSq(m23sq,m13sq))>0.01 ){
			std::cout<<m12sq<<" "<<kin->getThirdVariableSq(m23sq,m13sq)<<std::endl;
			std::cout<<"   " <<m23sq<<" "<<m13sq<<" "<<m12sq<<std::endl;
		}
		TParticle fparticleGam(22,1,0,0,0,0,*pGamma,W);
		TParticle fparticlePip(211,1,0,0,0,0,*pPip,W);
		TParticle fparticlePim(-211,1,0,0,0,0,*pPim,W);

		//call physics module
//		vector<double> x;
//		x.push_back(sqrt(m23sq));
//		x.push_back(sqrt(m13sq));
//		x.push_back(sqrt(m12sq));
		//ParameterList intensL = testBW.intensity(dataP, minPar);
		double AMPpdf = testBW.intensity(dataP).GetParameterValue(0);
		//double AMPpdf = testBW.intensity(x, minPar);

		double test = rando.Uniform(0,maxTest);
		double testmc = rando.Uniform(0,1.);

		//mb.setVal(m13);
		//double m13pdf = totAmp13.getVal();//fun_combi2->Eval(m13);
		if(maxTest<(weight*AMPpdf))
			cout << "Einschwingen zu kurz!" << endl;
		if(i<MaxEvents && test<(weight*AMPpdf)){
			if(outCnt==maxCnt){
				outCnt=0;
				cout << (i/(double)MaxEvents*100.) << "% : " << test << " " << (weight*AMPpdf) << endl;
			}
			outCnt++;
			i++;
			new((*fEvt)[0]) TParticle(fparticleGam);
			new((*fEvt)[1]) TParticle(fparticlePip);
			new((*fEvt)[2]) TParticle(fparticlePim);

			fTree.Fill();
		}

		if(mc<MaxEvents && testmc<weight){
			mc++;
			new((*fEvtPHSP)[0]) TParticle(fparticleGam);
			new((*fEvtPHSP)[1]) TParticle(fparticlePip);
			new((*fEvtPHSP)[2]) TParticle(fparticlePim);

			fTreePHSP.Fill();
		}
	}while(i<MaxEvents || mc<MaxEvents);
	cout << "100%! Write Data" << endl;

	fTree.Print();
	fTree.Write();
	fTreePHSP.Write();
	output.Close();

	testBW.printAmps();
	cout << "Done ... " << maxTest << endl << endl;

	return 0;
}
