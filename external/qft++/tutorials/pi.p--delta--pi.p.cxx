/* Copyright 2008 Mike Williams (mwill@jlab.org)
 *
 * This file is part of qft++.
 *
 * qft++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * qft++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with qft++.  If not, see <http://www.gnu.org/licenses/>.
 */
//_____________________________________________________________________________
/** @file pi.p--delta--pi.p.cxx
 *  @author Mike Williams
 *
 *  @brief Calculates intensity for \f$p\pi\rightarrow\Delta\rightarrow p\pi\f$
 *
 *  This tutorial calculates the scattering amplitude for 
 *  \f$p\pi\rightarrow\Delta\rightarrow p\pi\f$ using the qft++ classes.
 *
 *  The amplitude calculated is
 *  \f[\bar{u}(p_f,m_f)p_f^{\mu} 
 *  2m_{\Delta}P^{\frac{3}{2}}_{\mu\nu}(p_{\Delta}) p_i^{\nu}
 *  u(p_i,m_i) BW(p_{\Delta})\f]
 *
 */
//_____________________________________________________________________________

#include <getopt.h>
#include "relativistic-quantum-mechanics.h"

using namespace std;
//_____________________________________________________________________________

/// Prints program usage to the screen
void PrintUsage();

/// Prints a line accross the screen
void PrintLine(char __c){ 
  for(int i = 0; i < 80; i++) cout << __c; 
  cout << endl;
}
//_____________________________________________________________________________

int main(int __argc,char *__argv[]){

  /*__________________________Parse the Command Line_________________________*/
  int c;
  while((c = getopt(__argc,__argv,"h")) != -1){
    switch(c){
    case 'h': // help option
      PrintUsage();
      return EXIT_SUCCESS;
      break;
    default:
      break;
    }
  }  
  double w; // center-of-mass energy
  double theta; // proton theta_cm
  double PI = 4*atan(1.);
  double mp = 0.93827,mpi = 0.13957;
  double p,e;

  /*__________________________________Set Up_________________________________*/
  Vector4<double> pi,pf,ptot; // inital,final proton momenta
  DiracGamma gamma; // gamma^{mu}
  DiracSpinor ui,uf; // initial, final p spinors
  DiracSpinor delta(3/2.); // delta spinor
  Spin mz_i,mz_f; // z spin projections
  Matrix<complex<double> > amp(1,1); // amplitude

  /*______________________________Set 4-Momenta______________________________*/
  PrintLine(':');
  cout << "Enter center-of-mass energy (GeV): ";
  cin >> w;
  cout << "Enter final state theta_cm (degrees): ";
  cin >> theta;
  theta *= PI/180.;
  e = (w*w + mp*mp - mpi*mpi)/(2*w);
  p = sqrt(e*e - mp*mp);

  pi.SetP4(sqrt(mp*mp + p*p),0,0,p);
  pf.SetP4(sqrt(mp*mp + p*p),p*sin(theta),0,p*cos(theta));
  ptot.SetP4(w,0,0,0);

  /*_________________Initialize Spinors/Polarization Vectors_________________*/
  ui.SetP4(pi,mp);
  uf.SetP4(pf,mp);
  delta.SetP4(ptot,1.232);

  /*___________________________Calculate Intensity___________________________*/
  double intensity = 0.;
  for(mz_i = -1/2.; mz_i <= 1/2.; mz_i++){ // sum over initial p spins
    for(mz_f = -1/2.; mz_f <= 1/2.; mz_f++){ // sum over final p spins
      amp = Bar(uf(mz_f))*pf*delta.PropagatorBW(0.12)*pi*ui(mz_i);
      intensity += norm(amp(0,0));
    }
  }
  // average over initial spins
  intensity /= 2.;  

  cout << "Intensity: " << intensity << endl;

  PrintLine(':');
  return EXIT_SUCCESS;
}
//_____________________________________________________________________________

void PrintUsage(){ 
  cout << "Usage: pi.p--delta--pi.p " << endl;
  cout << "This executable calculates the scattering intensity for pi p --> "
       << "Delta --> pi p\nfor a given kinematics." << endl;
}
//_____________________________________________________________________________
