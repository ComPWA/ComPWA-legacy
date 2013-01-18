#include "Utils.h"
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
/** @file Utils.C
 *  @brief Utility function source file
 */
//_____________________________________________________________________________
/// Utility function used by Clebsch 
inline double dfact(double __x){
  if((__x < 0.00001) && (__x >= 0.0)) return 1.;
  if(__x < 0) return 0.;
  return __x*dfact(__x - 1.);
}
//_____________________________________________________________________________
/// Returns i!
inline int factorial(int __i) {
  int f = 1;
  if((__i == 0)||(__i == 1)) f = 1;  
  else{
    while(__i > 0){
      f = f*__i;
      __i--;
    }
  }
  return f;
} 
//_____________________________________________________________________________
// define a few macros
#define MAX(x,y) (x>y ? x : y)
#define MIN(x,y) (x<y ? x : y)
//_____________________________________________________________________________

vector<LS> GetValidLS(const Spin &__j,int __parity,const Spin &__s1,int __p1,
		      const Spin &__s2,int __p2){

  vector<LS> valid_ls;
  LS ls;

  for(Spin S = abs(__s1 - __s2); S <= (__s1 + __s2); S++){
    for(int L = (int)abs(__j - S); L <= (int)(__j + S); L++){
      // check the parity
      if(abs(__p1*__p2*pow(-1.,L) - __parity) < 1.e-5) {
	ls.L = L;
	ls.S = S;
	valid_ls.push_back(ls);
      }
    }
  }
  return valid_ls;
}
//_____________________________________________________________________________

double Gamma(double __z){
  
  float x = abs(__z);
  // Coefficients for the series expansion
  double c[7] = { 2.5066282746310005, 76.18009172947146, -86.50532032941677
	       ,24.01409824083091, -1.231739572450155, 0.1208650973866179e-2
               ,-0.5395239384953e-5};

   double y = x;
   double tmp = x + 5.5;
   tmp = (x+0.5)*log(tmp)-tmp;
   double ser = 1.000000000190015;
   for(int i = 1; i < 7; i++) {
      y += 1;
      ser += c[i]/y;
   }
   double lnGamma = tmp+log(c[0]*ser/x); // natural log of Gamma(z)
   if(__z > 0) return exp(lnGamma);   
   else if(__z == 0.) return 0.;   
   else {
     double pi = 3.1415926535;     
     return (-1.*pi)/(exp(lnGamma)*sin(pi*x)*x);
   }   
}
//_____________________________________________________________________________

double Clebsch(const Spin &__j1,const Spin &__m1,const Spin &__j2,
	       const Spin &__m2,const Spin &__J,const Spin &__M){

  if((__m1 + __m2) != __M) return 0.;
  if(abs(__m1) > __j1) return 0;
  if(abs(__m2) > __j2) return 0;

  // convert to pure integers (each 2*spin)
  int j1 = (int)(2.*__j1);
  int m1 = (int)(2.*__m1);
  int j2 = (int)(2.*__j2);
  int m2 = (int)(2.*__m2);
  int J = (int)(2.*__J);
  int M = (int)(2.*__M);

  double n0,n1,n2,n3,n4,n5,d0,d1,d2,d3,d4,A,exp;
  int nu = 0;
  
  double sum = 0;
  while(((d3=(j1-j2-M)/2+nu) < 0)||((n2=(j1-m1)/2+nu) < 0 )) { nu++;}
  while (((d1=(J-j1+j2)/2-nu) >= 0) && ((d2=(J+M)/2-nu) >= 0) 
	 &&((n1=(j2+J+m1)/2-nu) >= 0 )){
    d3=((j1-j2-M)/2+nu);
    n2=((j1-m1)/2+nu);
    d0=dfact((double) nu);
    exp=nu+(j2+m2)/2;
    n0 = (double) pow(-1.,exp);
    sum += ((n0*dfact(n1)*dfact(n2))/(d0*dfact(d1)*dfact(d2)*dfact(d3)));
    nu++;
  }

  if (sum == 0) return 0;

  n0 = J+1;
  n1 = dfact((double) (J+j1-j2)/2);
  n2 = dfact((double) (J-j1+j2)/2);
  n3 = dfact((double) (j1+j2-J)/2);
  n4 = dfact((double) (J+M)/2);
  n5 = dfact((J-M)/2);
  
  d0 = dfact((double) (j1+j2+J)/2+1);
  d1 = dfact((double) (j1-m1)/2);
  d2 = dfact((double) (j1+m1)/2);
  d3 = dfact((double) (j2-m2)/2);
  d4 = dfact((double) (j2+m2)/2);
  
  A = ((double) (n0*n1*n2*n3*n4*n5))/((double) (d0*d1*d2*d3*d4));


  return sqrt(A)*sum;		
}
//_____________________________________________________________________________

double Wigner_d(const Spin &__j,const Spin &__m,const Spin &__n,double __beta){

  int J = (int)(2.*__j);
  int M = (int)(2.*__m);
  int N = (int)(2.*__n);
  int temp_M, k, k_low, k_hi;
  double const_term = 0.0, sum_term = 0.0, d = 1.0;
  int m_p_n, j_p_m, j_p_n, j_m_m, j_m_n;
  int kmn1, kmn2, jmnk, jmk, jnk;
  double kk;

  if (J < 0 || abs (M) > J || abs (N) > J) {
    cerr << endl;
    cerr << "d: you have entered an illegal number for J, M, N." << endl;
    cerr << "Must follow these rules: J >= 0, abs(M) <= J, and abs(N) <= J." 
	 << endl;
    cerr << "J = " << J <<  " M = " << M <<  " N = " << N << endl;
    return 0.;
  }
  
  if (__beta < 0) {
    __beta = fabs (__beta);
    temp_M = M;
    M = N;
    N = temp_M;
  }

  m_p_n = (M + N) / 2;
  j_p_m = (J + M) / 2;
  j_m_m = (J - M) / 2;
  j_p_n = (J + N) / 2;
  j_m_n = (J - N) / 2;
  
  kk = (double)factorial(j_p_m)*(double)factorial(j_m_m)
    *(double)factorial(j_p_n) * (double)factorial(j_m_n) ;
  const_term = pow((-1.0),(j_p_m)) * sqrt(kk);	
  
  k_low = MAX(0, m_p_n);
  k_hi = MIN(j_p_m, j_p_n);

  for (k = k_low; k <= k_hi; k++) {
    
    kmn1 = 2 * k - (M + N) / 2;
    jmnk = J + (M + N) / 2 - 2 * k;
    jmk = (J + M) / 2 - k;
    jnk = (J + N) / 2 - k;
    kmn2 = k - (M + N) / 2;
	
    sum_term += pow ((-1.0), (k)) *
      ((pow (cos (__beta / 2.0), kmn1)) * (pow (sin (__beta / 2.0), jmnk))) /
      (factorial (k) * factorial (jmk) * factorial (jnk) * factorial (kmn2));
  }

  d = const_term * sum_term;
  return d;
}
//_____________________________________________________________________________

complex<double> ReggePropagator(double __t,double __s,double __a,double __b,
				const Spin &__spin,int __sig,int __exp_fact){
 
  complex<double> prop,numerator,denominator;
  double alpha = __a*__t + __b; // trajectory
  double pi = 3.1415926535; 
  complex<double> i(0.,1.); 
  numerator = pow(__s,alpha - __spin)*pi*__a;
  numerator *= ((double)__sig + ((double)__exp_fact)*exp(-i*pi*alpha));
  double gamma_arg = alpha + 1. - __spin;
  if(gamma_arg < 0.) denominator = 2*pi/(abs(gamma_arg)*Gamma(abs(gamma_arg)));
  else if(gamma_arg == 0.) denominator = 2*pi;
  else denominator = 2*sin(pi*gamma_arg)*Gamma(gamma_arg);

  prop = numerator/denominator;

  return prop;
}
//_____________________________________________________________________________

Spin GetSpin(const string &__spin){
  
  Spin j;
  string::size_type div_pos = __spin.find("/");
  if(div_pos != string::npos){
    string tmp;
    tmp.assign(__spin,0,div_pos);
    double numer = atof(tmp.c_str());
    tmp.assign(__spin,div_pos + 1,__spin.size() - div_pos - 1);
    double denom = atof(tmp.c_str());
    j = numer/denom;
  }
  else j = atof(__spin.c_str());
  
  return j;
}
//_____________________________________________________________________________

void Wigner_d(const Spin &__jmax,double __beta,
	      map<Spin,map<Spin,map<Spin,double> > > &__d){
  __d.clear();
  Spin jmin,one_half = 1/2.;
  if(__jmax == 0.){
    __d[0][0][0] = 1.;
    return;
  }
  double cb = cos(__beta);
  double sb = sin(__beta);
  // j=1 d's
  map<Spin,map<Spin,double> > d1;
  d1[1][1] = (1+cb)/2.;
  d1[1][0] = -sb/sqrt(2.);
  d1[1][-1] = (1-cb)/2.;
  d1[0][1] = sb/sqrt(2.);
  d1[0][0] = cb;
  d1[0][-1] = -sb/sqrt(2.);
  d1[-1][1] = (1-cb)/2.;
  d1[-1][0] = sb/sqrt(2.);
  d1[-1][-1] = (1+cb)/2.;

  if(__jmax.Denominator() == 1){ // integral spins
    __d[0][0][0] = 1.0;
    if(__jmax == 0.) return;
    __d[1][1][1] = d1[1][1];
    __d[1][1][0] = d1[1][0];
    __d[1][1][-1] = d1[1][-1];
    __d[1][0][1] = d1[0][1];
    __d[1][0][0] = d1[0][0];
    __d[1][0][-1] = d1[0][-1];
    __d[1][-1][1] = d1[-1][1];
    __d[1][-1][0] = d1[-1][0];
    __d[1][-1][-1] = d1[-1][-1];
    if(__jmax == 1.) return;
    jmin = 2.;
  }
  else { // half-integral spins      
    __d[one_half][one_half][one_half] = cos(__beta/2.);
    __d[one_half][one_half][-one_half] = -sin(__beta/2.);      
    __d[one_half][-one_half][one_half] = sin(__beta/2.);
    __d[one_half][-one_half][-one_half] = cos(__beta/2.);
    if(__jmax == one_half) return;
    jmin = 3/2.;
  }

  for(Spin j = jmin; j <= __jmax; j++){
    for(Spin m = -j; m <= j; m++){
      for(Spin n = -j; n <= j; n++){
	double djmn = 0.;
	for(Spin mm = -1; mm <= 1; mm++){
	  for(Spin nn = -1; nn <= 1; nn++){
	    djmn += Clebsch(j-1,m-mm,1,mm,j,m)*Clebsch(j-1,n-nn,1,nn,j,n)
	      *__d[j-1][m-mm][n-nn]*d1[mm][nn];
	  }
	}
	__d[j][m][n] = djmn;
      }
    }
  }  
}
//_____________________________________________________________________________

void Wigner_D(const Spin &__jmax,double __alpha,double __beta,double __gamma,
	      map<Spin,map<Spin,map<Spin,complex<double> > > > &__D){

  complex<double> i(0.,1.);
  map<Spin,map<Spin,map<Spin,double> > > d;
  Wigner_d(__jmax,__beta,d);
  Spin jmin;
  if(d.find(0) != d.end()) jmin = 0;
  else jmin = 1/2.;
  for(Spin j = jmin; j <= __jmax; j++){
    for(Spin m = -j; m <= j; m++){
      for(Spin n = -j; n <= j; n++) 
	__D[j][m][n] = exp(-i*(m*__alpha + n*__gamma))*d[j][m][n];
    }
  }
}
//_____________________________________________________________________________
