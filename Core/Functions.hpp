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
#include <complex>
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

  //! Pure Virtual interface for executing a strategy
  virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter> out) = 0;
};

class AddAll : public Strategy
{
public:
  AddAll(){
  };

  virtual const std::string to_str() const{
    return "+";
  }

  virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter> out){
    out = std::shared_ptr<AbsParameter>();
    unsigned int nMC = paras.GetNMultiComplex();
    unsigned int nMD = paras.GetNMultiDouble();
    unsigned int nC = paras.GetNComplex();
    unsigned int nD = paras.GetNDouble();
    unsigned int nI = paras.GetNInteger();

    if(nMC+nMD+nD+nI==0){
      //TODO: exception no input
      return false;
    }

	switch(out->type()){

      case ParType::MCOMPLEX:{
        //output multi complex: treat everything non-complex as real, there must be multi complex input
        if(!nMC){
          //TODO: exception wrong input
          return false;
        }

        unsigned int nElements = paras.GetMultiComplex(0)->GetNValues();

        std::complex<double> result(0,0);//sum up all 1-dim input
        //sum up complex parameter
        for(unsigned int i=0; i<nC; i++){
          result+=paras.GetComplexParameter(i)->GetValue();
        }
        //sum up double parameter
        for(unsigned int i=0; i<nD; i++){
          result+=paras.GetDoubleParameter(i)->GetValue();
        }
        //sum up integer parameter
        for(unsigned int i=0; i<nI; i++){
          result+=paras.GetIntegerParameter(i)->GetValue();
        }

        //fill MultiComplex parameter
        std::vector<std::complex<double>> results(nElements, result);
        for(unsigned int i=0; i<nMD; i++){
          std::shared_ptr<MultiDouble> tmp = paras.GetMultiDouble(i);
          for(unsigned int ele=0; ele<tmp->GetNValues(); ele++)
            results[ele]+=tmp->GetValue(ele);
        }
        for(unsigned int i=0; i<nMC; i++){
          std::shared_ptr<MultiComplex> tmp = paras.GetMultiComplex(i);
          for(unsigned int ele=0; ele<tmp->GetNValues(); ele++)
            results[ele]+=tmp->GetValue(ele);
        }

        out = std::shared_ptr<AbsParameter>(new MultiComplex("AddAllResult",results));

        break;
      }//end multi double

	  case ParType::MDOUBLE:{
	    //output multi double: ignore complex pars, there must be multi double input
	    if(!nMD){
	      //TODO: exception wrong input
	      return false;
	    }
		unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();
		double result=0;//sum up all 1-dim input
        //sum up double parameter
        for(unsigned int i=0; i<nD; i++){
          result+=paras.GetDoubleParameter(i)->GetValue();
        }
        //sum up integer parameter
        for(unsigned int i=0; i<nI; i++){
          result+=paras.GetIntegerParameter(i)->GetValue();
        }
        //fill MultiDouble parameter
        std::vector<double> results(nElements, result);
        for(unsigned int i=0; i<nMD; i++){
          std::shared_ptr<MultiDouble> tmp = paras.GetMultiDouble(i);
          for(unsigned int ele=0; ele<tmp->GetNValues(); ele++)
            results[ele]+=tmp->GetValue(ele);
        }
		out = std::shared_ptr<AbsParameter>(new MultiDouble("AddAllResult",results));

		break;
	  }//end multi double

      case ParType::COMPLEX:{
        //output complex: collapse everything non-complex as real-part
        std::complex<double> result(0,0);

        //sum up complex parameter
        for(unsigned int i=0; i<nC; i++){
          result+=paras.GetComplexParameter(i)->GetValue();
        }
        //sum up double parameter
        for(unsigned int i=0; i<nD; i++){
          result+=paras.GetDoubleParameter(i)->GetValue();
        }
        //sum up integer parameter
        for(unsigned int i=0; i<nI; i++){
          result+=paras.GetIntegerParameter(i)->GetValue();
        }
        //collapse MultiComplex parameter
        for(unsigned int i=0; i<nMC; i++){
          std::shared_ptr<MultiComplex> tmp = paras.GetMultiDouble(i);
          for(unsigned int ele=0; ele<tmp->GetNValues(); ele++)
            result+=tmp->GetValue(ele);
        }
        //collapse MultiDoubles parameter
        for(unsigned int i=0; i<nMD; i++){
          std::shared_ptr<MultiDouble> tmp = paras.GetMultiDouble(i);
          for(unsigned int ele=0; ele<tmp->GetNValues(); ele++)
            result+=tmp->GetValue(ele);
        }

        out = std::shared_ptr<AbsParameter>(new DoubleParameter("AddAllResult",result));
        break;
      }//end double

	  case ParType::DOUBLE:{
	    //output double: ignore complex pars, collapse everything else
	    double result=0;

	    //sum up double parameter
	    for(unsigned int i=0; i<nD; i++){
	      result+=paras.GetDoubleParameter(i)->GetValue();
	    }
        //sum up integer parameter
        for(unsigned int i=0; i<nI; i++){
          result+=paras.GetIntegerParameter(i)->GetValue();
        }
        //collapse MultiDoubles parameter
        for(unsigned int i=0; i<nMD; i++){
          std::shared_ptr<MultiDouble> tmp = paras.GetMultiDouble(i);
          for(unsigned int ele=0; ele<tmp->GetNValues(); ele++)
            result+=tmp->GetValue(ele);
        }

	    out = std::shared_ptr<AbsParameter>(new DoubleParameter("AddAllResult",result));
	    break;
	  }//end double

	  default:{
	    //TODO: exception output partype wrong
	    return false;
	  }

	}//end switch

	//return std::shared_ptr<AbsParameter>(new DoubleParameter("MultAllResult",result));
	return true;
  };
};

class MultAll : public Strategy
{
public:
  MultAll(){
  };

  virtual const std::string to_str() const{
    return "*";
  };

  virtual bool execute(ParameterList& paras,  std::shared_ptr<AbsParameter> out){
    out = std::shared_ptr<AbsParameter>();
    unsigned int nMC = paras.GetNMultiComplex();
    unsigned int nMD = paras.GetNMultiDouble();
    unsigned int nD = paras.GetNDouble();
    unsigned int nI = paras.GetNInteger();

    if(nMC+nMD+nD+nI==0){
      //TODO: exception no input
      return false;
    }

    switch(out->type()){

      case ParType::MDOUBLE:{
        //output multi double: ignore complex pars, there must be multi double input
        if(!nMD){
          //TODO: exception wrong input
          return false;
        }
        unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();
        double result=1.;//sum up all 1-dim input
        //mult up double parameter
        for(unsigned int i=0; i<nD; i++){
          result*=paras.GetDoubleParameter(i)->GetValue();
        }
        //mult up integer parameter
        for(unsigned int i=0; i<nI; i++){
          result*=paras.GetIntegerParameter(i)->GetValue();
        }
        //fill MultiDouble parameter
        std::vector<double> results(nElements, result);
        for(unsigned int i=0; i<nMD; i++){
          std::shared_ptr<MultiDouble> tmp = paras.GetMultiDouble(i);
          for(unsigned int ele=0; ele<tmp->GetNValues(); ele++)
            results[ele]*=tmp->GetValue(ele);
        }
        out = std::shared_ptr<AbsParameter>(new MultiDouble("MultAllResult",results));

        break;
      }//end multi double

      case ParType::DOUBLE:{
        //output double: ignore complex pars, collapse everything else
        double result=1.;

        //mult up double parameter
        for(unsigned int i=0; i<nD; i++){
          result*=paras.GetDoubleParameter(i)->GetValue();
        }
        //mult up integer parameter
        for(unsigned int i=0; i<nI; i++){
          result*=paras.GetIntegerParameter(i)->GetValue();
        }
        //collapse MultiDoubles parameter
        for(unsigned int i=0; i<nMD; i++){
          std::shared_ptr<MultiDouble> tmp = paras.GetMultiDouble(i);
          for(unsigned int ele=0; ele<tmp->GetNValues(); ele++)
            result*=tmp->GetValue(ele);
        }

        out = std::shared_ptr<AbsParameter>(new DoubleParameter("MultAllResult",result));
        break;
      }//end double

      default:{
        //TODO: exception output partype wrong
        return false;
      }

    }//end switch

    return true;
  };
};

class LogOf : public Strategy
{
public:
  LogOf(){
  };

  virtual const std::string to_str() const{
    return "Log";
  };

  virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter> out){
    out = std::shared_ptr<AbsParameter>();
    unsigned int nMC = paras.GetNMultiComplex();
    unsigned int nMD = paras.GetNMultiDouble();
    unsigned int nD = paras.GetNDouble();
    unsigned int nI = paras.GetNInteger();

    if(nMC+nMD+nD+nI==0){
      //TODO: exception no input
      return false;
    }
    //only one parameter possible
    if( (nMC+nMD+nD+nI)>1 ){
      //TODO: exception wrong input
      return false;
    }

    switch(out->type()){

      case ParType::MDOUBLE:{
        //output multi double: input must be one multi double
        if(!nMD){
          //TODO: exception wrong input
          return false;
        }
        unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();
        //fill MultiDouble parameter
        std::vector<double> results(nElements, 0.);
        std::shared_ptr<MultiDouble> tmp = paras.GetMultiDouble(0);
        for(unsigned int ele=0; ele<tmp->GetNValues() || ele<nElements; ele++)
          results[ele] = log(tmp->GetValue(ele));

        out = std::shared_ptr<AbsParameter>(new MultiDouble("LogResult",results));

        break;
      }//end multi double

      case ParType::DOUBLE:{
        //output double: log of one double input
        out = std::shared_ptr<AbsParameter>(new DoubleParameter("LogResult",log(paras.GetParameterValue(0))));
        break;
      }//end double

      default:{
        //TODO: exception output partype wrong
        return false;
      }

    }//end switch

    return true;
  };
};

class Square : public Strategy
{
public:
  Square(){
  };

  virtual const std::string to_str() const{
    return "^2";
  };

  virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter> out){
    out = std::shared_ptr<AbsParameter>();
    unsigned int nMC = paras.GetNMultiComplex();
    unsigned int nMD = paras.GetNMultiDouble();
    unsigned int nD = paras.GetNDouble();
    unsigned int nI = paras.GetNInteger();

    if(nMC+nMD+nD+nI==0){
      //TODO: exception no input
      return false;
    }
    //only one parameter possible
    if( (nMC+nMD+nD+nI)>1 ){
      //TODO: exception wrong input
      return false;
    }

    switch(out->type()){

      case ParType::MDOUBLE:{
        //output multi double: input must be one multi double
        if(!nMD){
          //TODO: exception wrong input
          return false;
        }
        unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();
        //fill MultiDouble parameter
        std::vector<double> results(nElements, 0.);
        std::shared_ptr<MultiDouble> tmp = paras.GetMultiDouble(0);
        for(unsigned int ele=0; ele<tmp->GetNValues() || ele<nElements; ele++)
          results[ele] = (tmp->GetValue(ele))*(tmp->GetValue(ele));

        out = std::shared_ptr<AbsParameter>(new MultiDouble("SqResult",results));

        break;
      }//end multi double

      case ParType::DOUBLE:{
        //output double: log of one double input
        out = std::shared_ptr<AbsParameter>(new DoubleParameter("SqResult",(paras.GetParameterValue(0)*paras.GetParameterValue(0))));
        break;
      }//end double

      default:{
        //TODO: exception output partype wrong
        return false;
      }

    }//end switch

    return true;
  };
};

class Power : public Strategy
{
public:
  Power(){
  };

  virtual const std::string to_str() const{
    return "^";
  }

  virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter> out){
    out = std::shared_ptr<AbsParameter>();
    unsigned int nMC = paras.GetNMultiComplex();
    unsigned int nMD = paras.GetNMultiDouble();
    unsigned int nD = paras.GetNDouble();
    unsigned int nI = paras.GetNInteger();

    if(nMC+nMD+nD+nI==0){
      //TODO: exception no input
      return false;
    }
    //only two double parameter possible
    if( !(nD+nMD==2) ){
      //TODO: exception wrong input
      return false;
    }

    switch(out->type()){

      case ParType::MDOUBLE:{
        //output multi double: input must be two multi double
        if(!(nMD==2)){
          //TODO: exception wrong input
          return false;
        }
        unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();
        //fill MultiDouble parameter
        std::vector<double> results(nElements, 0.);
        std::shared_ptr<MultiDouble> tmpA = paras.GetMultiDouble(0);
        std::shared_ptr<MultiDouble> tmpB = paras.GetMultiDouble(1);
        for(unsigned int ele=0; ele<nElements; ele++)
          results[ele] = pow(tmpA->GetValue(ele),tmpB->GetValue(ele));

        out = std::shared_ptr<AbsParameter>(new MultiDouble("PowResult",results));

        break;
      }//end multi double

      case ParType::DOUBLE:{
        //output double: power of two double input
        out = std::shared_ptr<AbsParameter>(new DoubleParameter("PowResult",pow(paras.GetParameterValue(0),paras.GetParameterValue(1))));
        break;
      }//end double

      default:{
        //TODO: exception output partype wrong
        return false;
      }

    }//end switch

    return true;
  };
};

#endif
