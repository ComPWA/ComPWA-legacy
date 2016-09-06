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
//! Plotting Tool for Physics-IF amplitudes.
/*!
 * @file AmpPlotterApp.cpp
 * A small GUI to plot amplitudes provided by Physics-Modules. Several models can
 * be illustrated with automatically generated sliders for their parameters.
*/

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

// Physics Interface header files go here
#include "Physics/BreitWigner/BreitWigner.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"

// Application header files
#include "test/TTripleSliderDemo.hpp"

// Root header files
#include <TApplication.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include "TGButton.h"
#include "TRootEmbeddedCanvas.h"
#include "TGLayout.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGTextEntry.h"
#include "TGSlider.h"
#include "TStyle.h"
#include "TGLabel.h"
#include "TGDoubleSlider.h"

#include "TComplex.h"
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TArrow.h"

using namespace std;

using ComPWA::Amplitude;
//using ComPWA::Physics::BreitWigner::BreitWigner;
using ComPWA::ParameterList;

using namespace ComPWA;

//namespace AmpPlotter{

  //TGMainFrame *fMain;

  enum ETestCommandIdentifiers {
     HId1,
     HId2,
     HId3,
     HCId1,
     HCId2,

     HSId1
  };


//}//end namespace AmpPlotter

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  std::cout << "  ComPWA Copyright (C) 2013  Mathias Michel " << std::endl;
  std::cout << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt" << std::endl;
  std::cout << std::endl;

  TApplication theApp("App", &argc, argv);
 // fMain = new TGMainFrame(gClient->GetRoot(),200,200);

  unsigned int whichAmp=0;
  shared_ptr<Amplitude> plotMe;
  //ParameterList par;
  unsigned int nPar;
  bool rel=false;
  double par[100];
  double slmin[100];
  double slmax[100];
  double step[100];
  TString name[100];

  switch(whichAmp){
    case 0:{ //Klaus 2BW-Summe
      nPar=7;
      double para[7] = {500, 2500,  150,  250, 2200,   50,   0};
      double slmina[7]= {   0,  300,    1,    0,  300,    1,   0};
      double slmaxa[7]= {1000, 3000,  300, 1000, 3000,  300, 360};
      double stepa[7]= {  10,  100,   10,   10,  100,   10,   2};
      name[0]= "A1";
      name[1]= "m1";
      name[2]= "G1";
      name[3]= "A2";
      name[4]= "m2";
      name[5]= "G2";
      name[6]= "Phase";
      rel=true;
      memcpy(par, para, nPar*sizeof(double));
      memcpy(slmin, slmina, nPar*sizeof(double));
      memcpy(slmax, slmaxa, nPar*sizeof(double));
      memcpy(step, stepa, nPar*sizeof(double));
      break;
    }
    case 1:{ //ComPWA single BW
      plotMe.reset(new ComPWA::Physics::BreitWigner::BreitWigner(0.,5.));
      ParameterList partmp;
      plotMe->copyParameterList(partmp);
      nPar = 0;
      for(unsigned int i=0; i<partmp.GetNDouble(); i++){
        std::shared_ptr<DoubleParameter> tmpdouble = partmp.GetDoubleParameter(i);
        if(!tmpdouble->IsFixed()){
          par[i] = tmpdouble->GetValue();
          slmin[i] = tmpdouble->GetMinValue();
          slmax[i] = tmpdouble->GetMaxValue();
          step[i] = 50;
          nPar++;
        }
      }
      break;
    }
  }

  //int ival0 = 500, ival1 = 2500, ival2 = 150, ival3 = 250, ival4 = 2200, ival5 = 50, ival6 =0; bool ival7=1;
  TTripleSliderDemo *plotGUI = new TTripleSliderDemo(par, slmin, slmax, step, name, nPar, rel);

  theApp.Run();

  cout << "Done ..." << endl << endl;

  return 0;
}
