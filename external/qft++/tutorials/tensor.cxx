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
/** @file tensor.cxx
 *  @author Mike Williams
 *
 *  @brief Provides example usage of the tensor module
 *
 *  This tutorial gives a number of examples of how to perform common tensor
 *  operations using the qft++ tensor package.
 */
//_____________________________________________________________________________

#include <getopt.h>
#include "tensor.h"

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
  //extern char* optarg;
  //extern int optind;
  
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

  /*_____________________________Creating Tensors____________________________*/
  // Tensor is a template class, here we'll work with double's and 
  // complex<double>'s...but you could use any data type that defines the 
  // proper operators. The constructor argument is the rank.
  Tensor<double> x1(1),x2(2),x3(3);
  Tensor<complex<double> > y1(1);
  Vector4<double> v;
  LeviCivitaTensor eps; // totally anti-symetric 4th rank tensor
  MetricTensor g; // g_{mu,nu} 

  /*_____________________________Setting Tensors_____________________________*/
  // The metric and Levi-Civita tensors are set when they're created. To set
  // elements of our other tensors we have several options. We can either set
  // them using other Tensor's (using the = operator), or by direct access to
  // the elements.

  PrintLine(':');
  cout << "We've set up the following tensors: " << endl;
  cout << "tensors storing double's:" << endl;
  // element access
  for(int e = 0; e < 4; e++) {
    x1(e) = e;
    x2(e,e) = 1.0; // make it diagnol
  }
  cout << "->x1: " << x1 << endl;
  cout << "->x2: " << x2 << endl;
 
  // using the assignment operator
  x3 = x1 % x2; // tensor outer product
  cout << "->x3: x1^{mu}x2^{nu,rho} (printing only defined for rank <= 2)" 
       << endl;
  // just make sure these don't throw errors
  x3.Symmetric();
  x3.AntiSymmetric();

  cout << "->metric: g^{mu,nu}: " << g << endl;
  cout << "->Levi-Civita: epsilon^{mu,nu,rho,sigma} (not printable)" 
       << endl;

  cout << "4-vectors storing double's:" << endl;
  double mass = 0.13957;
  v.SetP4(sqrt(1 + mass*mass),0,0,1.);
  cout << "->v: " << v << endl;

  cout << "tensors storing complex double's:" << endl;
  y1(0) = complex<double>(0.,1.);
  y1(3) = complex<double>(0.,1.);
  cout << "->y1: " << y1 << endl;

  /*_____________________________Basic Operations____________________________*/
  PrintLine(':');

  cout << "Tensor contractions:" << endl;

  // To contract the last index of x with the 1st index of y, just do x*y
  cout << "contraction of 1 index: " << endl;
  cout << "x1_{mu}x1^{mu}:\n->" << x1*x1 << endl;
  cout << "x1_{mu}x2^{mu,nu}:\n->" << x1 * x2 << endl;
  cout << "x1_{mu}x3^{mu,nu,rho}:\n->" << x1 * x3 << endl;
  cout << "x1_{rho}x3^{mu,nu,rho}:\n->" << x3 * x1 << endl;
  cout << "x1_{nu}x3^{mu,nu,rho}:\n->" << (x3.Permute(2,3)) * x1 << endl;
  cout << "sqrt(v^{mu}v_{mu}):\n->" << sqrt(v*v) << endl;

  // To contract the last n indicies of x with the 1st n indicies of y, do
  // x.Contract(y,n). If n is the rank of either x or y, then you can just do
  // (x|y)...a tensor inner product.
  cout << "contraction of multiple indicies:" << endl;
  cout << "x2_{mu,nu}x3^{mu,nu,rho}:\n->" << (x2|x3) << endl;
  cout << "x1^{mu}x2^{nu,rho}x3_{mu,nu,rho}:\n->" << ((x1%x2)|x3) << endl;
  cout << "(x1^{mu}x2^{nu,rho} + 2*x3^{mu,nu,rho})x3_{mu,nu,rho}:\n->"
       << (((x1%x2) + 2*x3)|x3) << endl;
  cout << "epsilon^{mu,nu,rho,pi} x2_{mu,nu}:\n->" << (eps|x2) << endl;
  cout << "etc...as complicated as you want, it's still easy in the code."
       << endl;

  PrintLine(':');  
  // We can also obtain symmeterized versions of tensors.
  cout << "Symmetrization is easy: " << endl;
  Tensor<double> xx(2);
  for(int i = 0; i < 4; i++) xx(i,i) = 1.0;
  xx(2,1) = -2.;
  xx(2,3) = 3.;
  cout << "for a tensor:\n->" << xx << endl;
  cout << "symmetric:\n->" << xx.Symmetric() << endl;
  cout << "anti-symmetric:\n->" << xx.AntiSymmetric() << endl;
  cout << "this works for any rank." << endl;

  // Lorentz transformations are done either thru a general transformation
  // tensor via x.Transform(lambda), where lambda_{mu,nu} is a 2nd rank
  // tensor...or to simply boost, you can use either a 1st rank tensor via
  // x.Boost(y) to boost to y's rest frame or specify the 3 boost vector
  // components and using x.Boost(bx,by,bz). You can also rotate by providing
  // the 3 euler angles via x.Rotate(alpha,beta,gamma).
  cout << "Lorentz transformations are also easy:" << endl;
  Vector4<double> v_cm(v);
  v_cm.Boost(v);
  cout << "boost to v's rest frame:\n->v_cm:" << v_cm << endl; 
  cout << "we can boost to any frame, rotate, etc..." << endl;

  /*_________________________Special 4-Vector Methods________________________*/
  PrintLine(':');
  cout << "There are also a number of special methods defined for 4-vectors:"
       << endl;
  cout << "4-vector v:->" << v << endl;
  cout << "invariant mass:-> " << v.Mass() << endl;
  cout << "beta:-> " << v.Beta() << endl;
  cout << "cos(theta):-> " << v.CosTheta() << endl;
  cout << "see the Vector4 class documentation for a complete list." << endl;
  PrintLine(':');

  return EXIT_SUCCESS;
}
//_____________________________________________________________________________

void PrintUsage(){ 
  cout << "Usage: tensor " << endl;
  cout << "This executable provides a number of example usages of the tensor "
       << "package. Run\nthe executable to see what's being done, then look at"
       << " the source file to see \nhow it's done in the code."
       << endl;
}
//_____________________________________________________________________________
