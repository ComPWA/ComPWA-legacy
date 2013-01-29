//! Test-Application of internal PWAParameter.
/*!
 * @file ParameterTestApp.cpp
 * This tiny application tests the ComPWA internal Parameter class PWAParameter.
*/

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

//Core header files go here
#include "Core/PWAParameter.hpp"
#include "Core/PWAParameterList.hpp"
#include "Core/PWAGenericPar.hpp"

using namespace std;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  //Test constructors
  PWAGenericPar<int> a; //empty
  PWAGenericPar<double> b; //empty double
  PWAGenericPar<int> c(2,0,5,1); //int par
  PWAGenericPar<int> d(c); //copy constructor
  PWAGenericPar<int> e(7,10,5,1); //contructor with wrong bounds
  shared_ptr<PWAGenericPar<int> > p(new PWAGenericPar<int>(3,0,5,1)); //pointer
  vector<PWAGenericPar<int> > v, w; //vector
  for(unsigned int par=0; par<10; par++)
    v.push_back(PWAGenericPar<int>(par,0,10,1));
  w = v; //copy vector

  //Test reading
  cout << "Constructor-Test output: \t" << endl;
  cout << "Empty int: \t" << a << endl;
  cout << "Empty double: \t" << b << endl;
  cout << "int: \t\t" << c << endl;
  cout << "copied int: \t" << d << endl;
  cout << "pointer int: \t" << *p << endl;
  for(unsigned int par=0; par<10; par++)
    cout << "vector " << par << ": \t" << v[par] << endl;
  for(unsigned int par=0; par<10; par++)
    cout << "copyvec " << par << ": \t" << w[par] << endl;
  cout << endl;

  //Test Getter & Setter
  cout << "Get & Set output: \t" << endl;
  cout << "Initial: \t" << a << endl;
  a.SetTValue(7); a.SetTMaxValue(10); a.SetTError(1);
  cout << "Setted: \t" << a << endl;
  cout << "GetVal: \t\t" << a.GetTValue() << endl;
  cout << "GetMin: \t\t" << a.GetTMinValue() << endl;
  cout << "GetMax: \t\t" << a.GetTMaxValue() << endl;
  cout << "GetErr: \t\t" << a.GetTError() << endl;
  cout << "Final: \t\t" << a << endl;
  cout << endl;

  //Test bound check
  cout << "bound output (e(7,10,5,1): \t" << endl;
  cout << "Initial: \t\t" << e << endl;
  cout << "HasBounds: \t\t\t" << e.HasBounds() << endl;
  e.SetTMaxValue(-1);
  cout << "After SetMaxValue(-1): \t" << e << endl;
  e.SetTMaxValue(6);
  cout << "After SetMaxValue(6): \t" << e << endl;
  e.SetTMinMax(10,5);
  cout << "After SetMinMax(10,5): \t" << e << endl;
  e.SetTMinMax(5,10);
  cout << "After SetMinMax(5,10): \t" << e << endl;
  cout << endl;

  //Test parameter list
  cout << endl;
  cout << "Test of PWAParameterList" << endl;
  PWAParameterList testempty;
  cout << "Empty Constructor: \t" << testempty << endl;
  PWAParameterList ints(v);
  cout << "IntVec Constructor: \t" << ints << endl;
  ints.AddParameter(b);
  cout << "IntVec added float: \t" << ints << endl;
  PWAGenericPar<int> toFill; ints.GetParameter(2,toFill);
  cout << "IntVec get int 2: \t" << toFill << endl;
  return 0;
}
