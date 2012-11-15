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
#include "PWAParameter.hpp"

using namespace std;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  //Test constructors
  PWAParameter<int> a; //empty
  PWAParameter<double> b; //empty double
  PWAParameter<int> c(2,0,5,1); //int par
  PWAParameter<int> d(c); //copy constructor
  PWAParameter<int> e(2,10,5,1); //contructor with wrong bounds
  shared_ptr<PWAParameter<int> > p(new PWAParameter<int>(3,0,5,1)); //pointer
  vector<PWAParameter<int> > v, w; //vector
  for(unsigned int par=0; par<10; par++)
    v.push_back(PWAParameter<int>(par,0,10,1));
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
    cout << "copyvec " << par << ": \t" << v[par] << endl;
  cout << endl;

  //Test Getter & Setter
  cout << "Get & Set output: \t" << endl;
  cout << "Initial: \t" << a << endl;
  a.SetValue(7); a.SetMaxValue(10); a.SetError(1);
  cout << "Setted: \t" << a << endl;
  cout << "GetVal: \t\t" << a.GetValue() << endl;
  cout << "GetMin: \t\t" << a.GetMinValue() << endl;
  cout << "GetMax: \t\t" << a.GetMaxValue() << endl;
  cout << "GetErr: \t\t" << a.GetError() << endl;
  cout << "Final: \t\t" << a << endl;
  cout << endl;

  //Test bound check
  cout << "bound output (e(2,10,5,1): \t" << endl;
  cout << "Initial: \t\t" << e << endl;
  cout << "HasBounds: \t\t\t" << e.HasBounds() << endl;
  e.SetMaxValue(-1);
  cout << "After SetMaxValue(-1): \t" << e << endl;
  e.SetMaxValue(6);
  cout << "After SetMaxValue(6): \t" << e << endl;
  e.SetMinMax(10,5);
  cout << "After SetMinMax(10,5): \t" << e << endl;
  e.SetMinMax(5,10);
  cout << "After SetMinMax(5,10): \t" << e << endl;
  cout << endl;

  return 0;
}
