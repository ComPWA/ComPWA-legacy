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
/** @file compton_scattering.cxx
 *  @author Mike Williams
 *
 *  @brief Calculates compton scattering \f$\frac{d\sigma}{dcos(\theta)}\f$.
 *
 *  This tutorial shows how to use qft++ classes to calculate the compton
 *  scattering differential cross section. Both the s- and u-channel 
 *  diagrams are calculated using the qft++ classes with no algebra needed.
 *  
 *  The amplitudes calculated are:
 *
 *  \f[A_{s-channel} = \bar{u}(p',m_p')(\gamma^{\mu}
 *  \epsilon^*_{\mu}(k',m_{\gamma}'))
 *  \frac{(p+k)^{\beta}\gamma_{\beta} + m}{(p+k)^2 -m^2}
 *  (\gamma^{\nu}\epsilon_{\nu}(k,m_{\gamma})) u(p,m_p)\f]
 *
 *  \f[A_{u-channel} = \bar{u}(p',m_p')(\gamma^{\nu}
 *  \epsilon_{\nu}(k,m_{\gamma}))
 *  \frac{(p-k')^{\beta}\gamma_{\beta} + m}{(p-k')^2 -m^2}
 *  (\gamma^{\mu}\epsilon^*_{\mu}(k',m_{\gamma}')) u(p,m_p)\f]
 *
 *  where p(p') denotes the initial(final) electron momentum and k(k') denotes
 *  the initial(final) photon momentum.
 *
 *  The user when running bin/compton_scattering is asked to provide the 
 *  kinematic info for the event, the incident photon energy along with the
 *  final photon angle. The qft++ answer is produced along with the famous
 *  Klein-Nishina equation result (which should give the same answer...if you
 *  find some set of kinematics where they don't email me). 
 *
 *  The user can then look at this source code to see how it calculates the
 *  amplitudes...once the 4-momenta are all set, it's 2 lines of code!
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
  double m_e = 0.511; // election mass
  double theta_f; // polar angle of final photon
  double e_gi,e_gf; // inital, final photon energy
  double PI = 4*atan(1.);
  double dsigma_dcos;

  /*__________________________________Set Up_________________________________*/
  Vector4<double> pi,pf,ki,kf; // inital,final e-,inital, final gamma
  DiracGamma gamma; // gamma^{mu}
  DiracSpinor ui,uf; // initial, final e- spinors
  DiracSpinor us,uu; // s-,u-channel exchange spinors
  PolVector epsi,epsf; // inital, final gamma polarization 4-vectors
  Spin mz_ei,mz_ef,mz_gi,mz_gf; // z spin projections
  Matrix<complex<double> > amp_s(1,1); // s-channel amplitude
  Matrix<complex<double> > amp_u(1,1); // u-channel amplitude
  complex<double> amp; // total amplitude

  /*______________________________Set 4-Momenta______________________________*/
  PrintLine(':');
  cout << "Enter incident photon energy (MeV): ";
  cin >> e_gi;
  cout << "Enter final photon angle (degrees): ";
  cin >> theta_f;
  theta_f *= PI/180.;
  e_gf = e_gi/(1 + e_gi*(1 - cos(theta_f))/m_e);

  ki.SetP4(e_gi,0,0,e_gi);
  pi.SetP4(m_e,0,0,0);
  kf.SetP4(e_gf,e_gf*sin(theta_f),0,e_gf*cos(theta_f));
  pf = pi + ki - kf;

  PrintLine(':');
  cout << "dsigma_dcos(theta) calculated by:" << endl;
  cout << "->Klein-Nishina equation: ";
  
  dsigma_dcos = (1/(16*m_e*m_e*PI))*pow(e_gf/e_gi,2)
    *(e_gf/e_gi + e_gi/e_gf - pow(sin(theta_f),2));
  cout << dsigma_dcos << endl;

  /*_________________Initialize Spinors/Polarization Vectors_________________*/
  ui.SetP4(pi,m_e);
  uf.SetP4(pf,m_e);
  epsi.SetP4(ki,0.);
  epsf.SetP4(kf,0.);
  us.SetP4(pi+ki,m_e);
  uu.SetP4(pi-kf,m_e);
  
  /*___________________________Calculate Intensity___________________________*/
  double intensity = 0.;
  for(mz_ei = -1/2.; mz_ei <= 1/2.; mz_ei++){ // sum over initial e- spins
    for(mz_gi = -1; mz_gi <= 1; mz_gi+=2){ // sum over initial gamma spins
      for(mz_ef = -1/2.; mz_ef <= 1/2.; mz_ef++){ // sum over final e- spins
	for(mz_gf = -1; mz_gf <= 1; mz_gf+=2){ // sum over final gamma spins
	  // s-channel piece
	  amp_s = Bar(uf(mz_ef))*(gamma*(epsf(mz_gf).Conjugate()))
	    *us.Propagator()*(gamma*epsi(mz_gi))*ui(mz_ei);
	  // u-channel piece
	  amp_u = Bar(uf(mz_ef))*(gamma*epsi(mz_gi))*uu.Propagator()
	    *(gamma*(epsf(mz_gf).Conjugate()))*ui(mz_ei);

	  amp = (amp_s + amp_u)(0,0); // total amp
	  intensity += norm(amp);
	}
      }
    }
  }
  // average over initial spins
  intensity /= 4.;  
  double dsigma_dcos_factor = (1./(2*m_e*2*e_gi))*(e_gf*e_gf/(8*PI*e_gi*m_e));
  dsigma_dcos = dsigma_dcos_factor*intensity;
  cout << "->qft++: " << dsigma_dcos << endl;

  PrintLine(':');
  return EXIT_SUCCESS;
}
//_____________________________________________________________________________

void PrintUsage(){ 
  cout << "Usage: compton_scattering " << endl;
  cout << "This executable calculates the compton scattering differential "
       << "cross section\nfor a given kinematics." << endl;
}
//_____________________________________________________________________________
