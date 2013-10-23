//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Functions to be used in FuntionTree.
/*! \class Strategy
 * \class AddAll
 * \class MultAll
 * \class PowerTwo
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
#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"
#include "Core/AbsParameter.hpp"

class Strategy
{
public:
  //! Constructor
  Strategy(){
  };

  //! friend function to stream parameter information to output
  /*!
   * Declaring the stream-operator << as friend allows to stream parameter
   * information to the output as easily as a generic type.
   * \sa make_str(), to_str()
  */
  friend std::ostream& operator<<( std::ostream& out, std::shared_ptr<Strategy> b ){
    return out << b->to_str();
  }

  //! friend function to stream parameter information to output
  /*!
   * Declaring the stream-operator << as friend allows to stream parameter
   * information to the output as easily as a generic type.
   * \sa make_str(), to_str()
  */
  friend std::ostream& operator<<( std::ostream& out, const Strategy& b ){
    return out << b.to_str();
  }

  //! Pure Virtual interface for streaming info about the strategy
  virtual const std::string to_str() const =0;

  //! Pure Virtual interface for executing a function
  virtual std::shared_ptr<AbsParameter> execute(const ParameterList& paras) = 0;
};

class AddAll : public Strategy
{
public:
  AddAll(){
  };

  virtual const std::string to_str() const{
    return "+";
  }

  virtual std::shared_ptr<AbsParameter> execute(const ParameterList& paras){
    double result = 0;
    for(unsigned int i=0; i<paras.GetNDouble(); i++)
      result+=paras.GetParameterValue(i);

    //ParameterList out;
    // out.AddParameter(DoubleParameter("AddAllResult",result));
    return std::shared_ptr<AbsParameter>(new DoubleParameter("AddAllResult",result));
  };
};

class MultAll : public Strategy
{
public:
  MultAll(){
  };

  virtual const std::string to_str() const{
    return "*";
  }

  virtual std::shared_ptr<AbsParameter> execute(const ParameterList& paras){
    double result = 1.;
    for(unsigned int i=0; i<paras.GetNDouble(); i++)
      result*=paras.GetParameterValue(i);

    //ParameterList out;
    //out.AddParameter(DoubleParameter("MultAllResult",result));
    return std::shared_ptr<AbsParameter>(new DoubleParameter("MultAllResult",result));
  };
};

class PowerTwo : public Strategy
{
public:
  PowerTwo(){
  };

  virtual const std::string to_str() const{
    return "^";
  }

  virtual std::shared_ptr<AbsParameter> execute(const ParameterList& paras){
    if(paras.GetNDouble()!=2){
        throw BadIndex("need exact two parameters");
        return NULL;
    }
    //return pow(paras[0],paras[1]);
    //ParameterList out;
    //out.AddParameter(DoubleParameter("PowerTwoResult",pow(paras.GetParameterValue(0),paras.GetParameterValue(1))));
    return std::shared_ptr<AbsParameter>(new DoubleParameter("PowerTwoResult",pow(paras.GetParameterValue(0),paras.GetParameterValue(1))));
  };
};

#endif
