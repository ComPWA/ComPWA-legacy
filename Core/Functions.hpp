//! Functions to be used in FuntionTree.
/*! \class Strategy
 * @file Functions.hpp
 * This file contains Functions implementing the Strategy interface so they
 * can be used inside a node of the FuntionTree to calculate the node-value.
 * In addition to the simple functions provided here, the interface can also
 * be used at other places to provide functions for the FunctionTree.
*/

#ifndef _FUNCTIONS_HPP_
#define _FUNCTIONS_HPP_

#include <vector>
#include <math.h>

#include "Core/Exceptions.hpp"

class Strategy
{
public:
  Strategy(){
  };

  virtual double execute(const std::vector<double>& paras) = 0;
};

class AddAll : public Strategy
{
public:
  AddAll(){
  };

  virtual double execute(const std::vector<double>& paras){
    double result = 0;
    for(unsigned int i=0; i<paras.size(); i++)
      result+=paras[i];
    return result;
  };
};

class MultAll : public Strategy
{
public:
  MultAll(){
  };

  virtual double execute(const std::vector<double>& paras){
    double result = 1.;
    for(unsigned int i=0; i<paras.size(); i++)
      result*=paras[i];
    return result;
  };
};

class PowerTwo : public Strategy
{
public:
  PowerTwo(){
  };

  virtual double execute(const std::vector<double>& paras){
    if(paras.size()!=2){
        throw BadIndex("need exact two parameters");
        return 0;
    }
    return pow(paras[0],paras[1]);
  };
};

#endif
