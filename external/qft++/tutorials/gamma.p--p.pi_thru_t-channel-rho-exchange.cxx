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
/** @file gamma.p--p.pi_thru_t-channel-rho-exchange.cxx
 *  @author Mike Williams
 *
 *  @brief Calculates intensity for \f$\gamma p\rightarrow p\pi\f$ with t-channel \f$\rho\f$ exchange.
 *
 *  This tutorial shows how to use the qft++ classes to calculate the 
 *  amplitude for \f$\gamma\pi\rightarrow p\pi\f$ where a \f$\rho\f$ meson is
 *  exchanged (t-channel).
 *
 *  The amplitude calculated is:
 *  \f[\bar{u}(p_f,m_f)\epsilon_{\mu\nu\rho\sigma}p^{\mu}_{\gamma}
 *   p^{\nu}_{\rho}\epsilon^{\rho}(p_{\gamma},m_{\gamma})(\gamma^{\sigma} +
 *   \frac{i\kappa_{\rho}}{2m_p}\sigma^{\sigma\alpha}p^{\rho}_{\alpha})
 *   u(p_i,m_i) BW(p_{\rho})\f]
 *
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
  double mp = 0.93827,mpi = 0.13957; // masses
  double theta; // polar angle of final pion
  double PI = 4*atan(1.);
  double w; // center-of-mass energy
  double p,e;
  double kappa_rho = 6.0; // guess a value...it doesn't matter for this
  complex<double> i(0,1.);

  /*__________________________________Set Up_________________________________*/
  Vector4<double> pi,pf,ppi,pg,prho; // 4-momenta
  DiracGamma gamma; // gamma^{mu}
  DiracSigma sigma; // sigma^{mu,nu}
  LeviCivitaTensor epsilon;
  DiracSpinor ui,uf; // initial, final p spinors
  PolVector eps_g; // polarization 4-vector for the gamma
  Spin mz_i,mz_f,mz_g; // z spin projections
  Matrix<complex<double> > amp(1,1); // amplitude

  /*______________________________Set 4-Momenta______________________________*/
  PrintLine(':');
  cout << "Enter center-of-mass energy (GeV): ";
  cin >> w;
  cout << "Enter final state theta_cm (degrees): ";
  cin >> theta;
  theta *= PI/180.;

  e = (w*w + mp*mp)/(2*w);
  p = sqrt(e*e - mp*mp);  
  pi.SetP4(sqrt(mp*mp + p*p),0,0,-p);
  pg.SetP4(p,0,0,p);

  e = (w*w + mp*mp - mpi*mpi)/(2*w);
  p = sqrt(e*e - mp*mp);  
  pf.SetP4(sqrt(mp*mp + p*p),-p*sin(theta),0,-p*cos(theta));
  ppi.SetP4(sqrt(mpi*mpi + p*p),p*sin(theta),0,p*cos(theta));

  prho = pg - ppi;

  /*_________________Initialize Spinors/Polarization Vectors_________________*/
  ui.SetP4(pi,mp);
  uf.SetP4(pf,mp);
  eps_g.SetP4(pg,0.);
  
  /*___________________________Calculate Intensity___________________________*/
  double intensity = 0.;
  for(mz_i = -1/2.; mz_i <= 1/2.; mz_i++){ // sum over initial p spins
    for(mz_g = -1; mz_g <= 1; mz_g+=2){ // sum over initial gamma spins
      for(mz_f = -1/2.; mz_f <= 1/2.; mz_f++){ // sum over final p spins

	amp = Bar(uf(mz_f))*(epsilon|(pg % prho % eps_g(mz_g))) 
	  * (gamma + (i*kappa_rho/(2*mp))*(sigma*prho)) * ui(mz_i);
	amp *= BreitWigner(prho,0.77,0.15); // add the rho's breit-wigner
	
	intensity += norm(amp(0,0));
      }
    }
  }

  // average over initial spins
  intensity /= 4.;  

  cout << "Intensity: " << intensity << endl;

  PrintLine(':');
  return EXIT_SUCCESS;
}
//_____________________________________________________________________________

void PrintUsage(){ 
  cout << "Usage: gamma.p--p.pi_thru_t-channel-rho-exchange " << endl;
  cout << "This executable calculates the scattering intensity for gamma p -->"
       << " pi p\nthru t-channel rho exchange for a given kinematics." << endl;
}
//_____________________________________________________________________________
