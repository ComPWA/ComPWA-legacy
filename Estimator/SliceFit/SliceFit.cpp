//-------------------------------------------------------------------------------
// Copyright (c) 2013 michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     michel - initial API and implementation
//-------------------------------------------------------------------------------
#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

#include "Estimator/SliceFit/SliceFit.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/ParameterList.hpp"

SliceFit::SliceFit(std::shared_ptr<Amplitude> inPIF, std::shared_ptr<Data> inDIF)
  : pPIF_(inPIF), pDIF_(inDIF)
   ,M(3.096916),Br(0.000093),m1(0),m2(0.139570),m3(0.139570),PI( 3.14159){

  Double_t m23_min = (m2+m3);
  Double_t m23_max = (M-m1);
  Double_t m13_min = (m1+m3);
  Double_t m13_max = (M-m2);
  Double_t m12_min = (m1+m2);
  Double_t m12_max = (M-m3);

  dalitzPlot_ = new TH2D("SliceDalitz","SliceDalitz", //TODO: variable binning
      100,m23_min*m23_min,m23_max*m23_max, 100,m13_min*m13_min,m13_max*m13_max);

  //fill dalitz
  double m23sq, m13sq, m12sq;

  for(unsigned int i=0; i<inDIF->GetNEvents();i++){
      Event event(myReader.getEvent(i));

      //if(!myReader.getEvent(i, event)) continue; TODO: try exception
      if(!event.getNParticles() == 3) continue;
      const Particle &a(event.getParticle(0));
      const Particle &b(event.getParticle(1));
      const Particle &c(event.getParticle(2));
      masssq12 = pow(a.E+b.E,2) - pow(a.px+b.px ,2) - pow(a.py+b.py ,2) - pow(a.pz+b.pz ,2);
      masssq13 = pow(a.E+c.E,2) - pow(a.px+c.px ,2) - pow(a.py+c.py ,2) - pow(a.pz+c.pz ,2);
      masssq23 = pow(b.E+c.E,2) - pow(b.px+c.px ,2) - pow(b.py+c.py ,2) - pow(b.pz+c.pz ,2);

    dalitzPlot_->Fill(masssq23 ,masssq13);
  }

}

SliceFit::~SliceFit(){

}

double SliceFit::controlParameter(ParameterList& minPar){
  unsigned int nEvents = pDIF_->getNEvents();
  unsigned int nParts = ((Event)pDIF_->getEvent(0)).getNParticles();

  //check if able to handle this many particles
  if(nParts!=3){
      //TODO: exception
      return 0;
  }

  //------------fit every slice

  return lh;
}
