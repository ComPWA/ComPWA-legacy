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
#include "Core/DataPoint.hpp"

class Strategy
{
public:
	//! Constructor
	Strategy(ParType in):checkType(in){
	};
	virtual ~Strategy(){}

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

	//! Get ParType
	virtual const ParType OutType() const {
		return checkType;
	}

	//! Pure Virtual interface for streaming info about the strategy
	virtual const std::string to_str() const =0;

	//! Pure Virtual interface for executing a strategy
	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) = 0;

protected:
	ParType checkType;
};

class Inverse : public Strategy
{
public:
	Inverse(ParType in):Strategy(in){
	};
	virtual ~Inverse() {};

	virtual const std::string to_str() const{
		return "+";
	}

	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out){
		if( checkType != out->type() ) return false;
		unsigned int nMC = paras.GetNMultiComplex();
		unsigned int nMD = paras.GetNMultiDouble();
		unsigned int nC = paras.GetNComplex();
		unsigned int nD = paras.GetNDouble();
		unsigned int nI = paras.GetNInteger();

		if(nMC+nC!=0){
			//TODO: can't handle complex numbers
			return false;
		}
		if(nMD+nD+nI!=1 ){
			//TODO: exception too many variables
			return false;
		}

		switch(checkType){
		case ParType::DOUBLE:{
			double var = paras.GetDoubleParameterValue(0);
			if(var==0){
				out = std::shared_ptr<AbsParameter>(new DoubleParameter( out->GetName(),0 ));
				//TODO: exception dividion by 0
				BOOST_LOG_TRIVIAL(error) << "Inverse strategy: division by 0";
			} else
				out = std::shared_ptr<AbsParameter>(new DoubleParameter( out->GetName(),1/var ));
			//    	  std::cout<<"Strategy Inverse (DOUBLE) "<<1/var<<std::endl;
			break;
		}//end double
		default:{
			//TODO: exception output partype wrong
			return false;
		}
		}//end switch
		return true;
	}//end execute
};

class SquareRoot : public Strategy
{
public:
	SquareRoot(ParType in):Strategy(in){
	};
	virtual ~SquareRoot(){}

	virtual const std::string to_str() const{
		return "+";
	}

	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out){
		if( checkType != out->type() ) return false;
		unsigned int nMC = paras.GetNMultiComplex();
		unsigned int nMD = paras.GetNMultiDouble();
		unsigned int nC = paras.GetNComplex();
		unsigned int nD = paras.GetNDouble();
		unsigned int nI = paras.GetNInteger();

		if(nMC+nC!=0){
			//TODO: can't handle complex numbers
			return false;
		}
		if(nMD+nD+nI!=1 ){
			//TODO: exception too many variables
			return false;
		}

		switch(checkType){
		case ParType::DOUBLE:{
			double var = paras.GetDoubleParameterValue(0);
			if(var<0){
				out = std::shared_ptr<AbsParameter>(new DoubleParameter( out->GetName(),-1 ));
				//TODO: exception argument <0
				BOOST_LOG_TRIVIAL(error) << "SquareRoot strategy: argument < 0! return -1";
			} else
				out = std::shared_ptr<AbsParameter>(new DoubleParameter( out->GetName(),sqrt(var) ));
			//    	  std::cout<<"Strategy SquareRoot (DOUBLE) "<<sqrt(var)<<std::endl;
			break;
		}//end double
		default:{
			//TODO: exception output partype wrong
			return false;
		}
		}//end switch
		return true;
	}//end execute
};
class AddAll : public Strategy
{
public:
	AddAll(ParType in):Strategy(in){
	};
	virtual ~AddAll(){}

	virtual const std::string to_str() const{
		return "+";
	}

	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out){
		if( checkType != out->type() ) return false;
		unsigned int nMC = paras.GetNMultiComplex();
		unsigned int nMD = paras.GetNMultiDouble();
		unsigned int nC = paras.GetNComplex();
		unsigned int nD = paras.GetNDouble();
		unsigned int nI = paras.GetNInteger();

		if(nMC+nMD+nD+nI+nC==0){
			//TODO: exception no input
			return false;
		}

		switch(checkType){

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

			out = std::shared_ptr<AbsParameter>(new MultiComplex(out->GetName(),results));

			break;
		}//end multi complex

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
			out = std::shared_ptr<AbsParameter>(new MultiDouble(out->GetName(),results));

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
				std::shared_ptr<MultiComplex> tmp = paras.GetMultiComplex(i);
				for(unsigned int ele=0; ele<tmp->GetNValues(); ele++)
					result+=tmp->GetValue(ele);
			}
			//collapse MultiDoubles parameter
			for(unsigned int i=0; i<nMD; i++){
				std::shared_ptr<MultiDouble> tmp = paras.GetMultiDouble(i);
				for(unsigned int ele=0; ele<tmp->GetNValues(); ele++)
					result+=tmp->GetValue(ele);
			}

			out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(),result));
			//        std::cout<<"Strategy AddAll (COMPLEX) "<<result<<std::endl;
			break;
		}//end complex

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

			//        std::cout<<"Strategy AddAll (DOUBLE) "<<result<<std::endl;
			out = std::shared_ptr<AbsParameter>(new DoubleParameter(out->GetName(),result));
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
	MultAll(ParType in):Strategy(in){
	};
	virtual ~MultAll(){}

	virtual const std::string to_str() const{
		return "*";
	};

	virtual bool execute(ParameterList& paras,  std::shared_ptr<AbsParameter>& out){
		if( checkType != out->type() ) return false;

		unsigned int nMC = paras.GetNMultiComplex();
		unsigned int nMD = paras.GetNMultiDouble();
		unsigned int nC = paras.GetNComplex();
		unsigned int nD = paras.GetNDouble();
		unsigned int nI = paras.GetNInteger();

		if(nMC+nMD+nD+nI+nC==0){
			//TODO: exception no input
			return false;
		}

		switch(checkType){

		case ParType::MCOMPLEX:{
			//output multi complex: treat everything non-complex as real, there must be multi complex input
			if(!nMC){
				//TODO: exception wrong input
				return false;
			}

			unsigned int nElements = paras.GetMultiComplex(0)->GetNValues();

			std::complex<double> result(1.,0.);//mult up all 1-dim input
			//mult up complex parameter
			for(unsigned int i=0; i<nC; i++){
				result*=paras.GetComplexParameter(i)->GetValue();
			}
			//mult up double parameter
			for(unsigned int i=0; i<nD; i++){
				result*=paras.GetDoubleParameter(i)->GetValue();
			}
			//mult up integer parameter
			for(unsigned int i=0; i<nI; i++){
				result*=paras.GetIntegerParameter(i)->GetValue();
			}

			//fill MultiComplex parameter
			std::vector<std::complex<double>> results(nElements, result);
			for(unsigned int i=0; i<nMD; i++){
				std::shared_ptr<MultiDouble> tmp = paras.GetMultiDouble(i);
				for(unsigned int ele=0; ele<tmp->GetNValues(); ele++)
					results[ele]*=tmp->GetValue(ele);
			}
			for(unsigned int i=0; i<nMC; i++){
				std::shared_ptr<MultiComplex> tmp = paras.GetMultiComplex(i);
				for(unsigned int ele=0; ele<tmp->GetNValues(); ele++)
					results[ele]*=tmp->GetValue(ele);
			}

			out = std::shared_ptr<AbsParameter>(new MultiComplex(out->GetName(),results));

			break;
		}//end multi complex

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
			out = std::shared_ptr<AbsParameter>(new MultiDouble(out->GetName(),results));

			break;
		}//end multi double

		case ParType::COMPLEX:{
			//output complex: collapse everything non-complex as real-part
			std::complex<double> result(1.,0);

			//mult up complex parameter
			for(unsigned int i=0; i<nC; i++){
				result*=paras.GetComplexParameter(i)->GetValue();
			}
			//mult up double parameter
			for(unsigned int i=0; i<nD; i++){
				result*=paras.GetDoubleParameter(i)->GetValue();
			}
			//mult up integer parameter
			for(unsigned int i=0; i<nI; i++){
				result*=paras.GetIntegerParameter(i)->GetValue();
			}
			//collapse MultiComplex parameter
			for(unsigned int i=0; i<nMC; i++){
				std::shared_ptr<MultiComplex> tmp = paras.GetMultiComplex(i);
				for(unsigned int ele=0; ele<tmp->GetNValues(); ele++)
					result*=tmp->GetValue(ele);
			}
			//collapse MultiDoubles parameter
			for(unsigned int i=0; i<nMD; i++){
				std::shared_ptr<MultiDouble> tmp = paras.GetMultiDouble(i);
				for(unsigned int ele=0; ele<tmp->GetNValues(); ele++)
					result*=tmp->GetValue(ele);
			}

			out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(),result));
			break;
		}//end complex

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

			out = std::shared_ptr<AbsParameter>(new DoubleParameter(out->GetName(),result));
			//          std::cout<<"Strategy MultAll (DOUBLE) "<<result<<std::endl;
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
	LogOf(ParType in):Strategy(in){
	};
	virtual ~LogOf(){};

	virtual const std::string to_str() const{
		return "Log";
	};

	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out){
		if( checkType != out->type() ) return false;

		//    unsigned int nMC = paras.GetNMultiComplex();
		unsigned int nMD = paras.GetNMultiDouble();
		//    unsigned int nC = paras.GetNComplex();
		unsigned int nD = paras.GetNDouble();
		//    unsigned int nI = paras.GetNInteger();

		if(nMD+nD==0){
			//TODO: exception no input
			return false;
		}
		//only one parameter possible
		if( (nMD+nD)>1 ){
			//TODO: exception wrong input
			return false;
		}

		switch(checkType){

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
				results[ele] = std::log(tmp->GetValue(ele));

			out = std::shared_ptr<AbsParameter>(new MultiDouble(out->GetName(),results));

			break;
		}//end multi double

		case ParType::DOUBLE:{
			if(!nD){
				//TODO: exception wrong input
				return false;
			}
			//output double: log of one double input
			out = std::shared_ptr<AbsParameter>(
					new DoubleParameter(
							out->GetName(),
							std::log(paras.GetDoubleParameterValue(0)
							)
					)
			);
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

class Complexify : public Strategy
{
public:
	Complexify(ParType in):Strategy(in){
	};
	virtual ~Complexify(){}

	virtual const std::string to_str() const{
		return "MakeComplex";
	};

	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out){
		if( checkType != out->type() ) return false;

		unsigned int nMC = paras.GetNMultiComplex();
		unsigned int nMD = paras.GetNMultiDouble();
		unsigned int nC = paras.GetNComplex();
		unsigned int nD = paras.GetNDouble();
		unsigned int nI = paras.GetNInteger();

		if(nMD+nD==0){
			//TODO: exception no input
			return false;
		}
		//only one parameter possible
		if( (nMD+nD)!=2 || (nMC+nC+nI)!=0 ){
			//TODO: exception wrong input
			return false;
		}

		switch(checkType){

		case ParType::MCOMPLEX:{
			//output multi complex: input must be two multi double
			if(!(nMD==2)){
				//TODO: exception wrong input
				return false;
			}
			unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();
			//fill MultiDouble parameter
			std::vector<std::complex<double> > results(nElements, std::complex<double>(0.,0.));
			std::shared_ptr<MultiDouble> ma = paras.GetMultiDouble(0);
			std::shared_ptr<MultiDouble> mphi = paras.GetMultiDouble(1);
			for(unsigned int ele=0; ele<nElements; ele++)
				results[ele] = std::complex<double>(std::fabs(ma->GetValue(ele))*std::cos(mphi->GetValue(ele)),
						std::fabs(ma->GetValue(ele))*std::sin(mphi->GetValue(ele)));//a*cos(phi),a*sin(phi)

			out = std::shared_ptr<AbsParameter>(new MultiComplex(out->GetName(),results));

			break;
		}//end multi complex

		case ParType::COMPLEX:{
			//output complex: input must be two double
			if(!(nD==2)){
				//TODO: exception wrong input
				return false;
			}
			double a = std::fabs(paras.GetDoubleParameter(0)->GetValue());
			double phi = paras.GetDoubleParameter(1)->GetValue();
			out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(),
					std::complex<double>(a*std::cos(phi),a*std::sin(phi))));
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

class AbsSquare : public Strategy
{
public:
	AbsSquare(ParType in):Strategy(in){
	};
	virtual ~AbsSquare(){}

	virtual const std::string to_str() const{
		return "||^2";
	};

	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out){
		if( checkType != out->type() ) return false;

		unsigned int nMC = paras.GetNMultiComplex();
		unsigned int nMD = paras.GetNMultiDouble();
		unsigned int nC = paras.GetNComplex();
		unsigned int nD = paras.GetNDouble();
		unsigned int nI = paras.GetNInteger();

		if(nMC+nMD+nD+nI+nC==0)
			throw std::runtime_error("AbsSquare Strategy: Input parameter list is empty!");

		//only one parameter possible
		if( (nMC+nMD+nD+nI+nC)>1 )
			throw std::runtime_error("AbsSquare Strategy: More then one type of parameter in input!");


		switch(checkType){
		case ParType::MDOUBLE:{
			//output multi double: input must be multi double or multi complex
			if(!nMD && !nMC)
				throw std::runtime_error("AbsSquare Strategy: Requested output is MDOUBLE but no MDOUBLE or MCOMPLEX was given as input!");
			if(nMD){
				unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();
				//fill MultiDouble parameter
				std::vector<double> results(nElements, 0.);
				std::shared_ptr<MultiDouble> tmp = paras.GetMultiDouble(0);
				for(unsigned int ele=0; ele<nElements; ele++)
					results[ele] = std::norm(tmp->GetValue(ele));

				out = std::shared_ptr<AbsParameter>(new MultiDouble(out->GetName(),results));
			}else if(nMC){
				unsigned int nElements = paras.GetMultiComplex(0)->GetNValues();
				//fill MultiDouble parameter
				std::vector<double> results(nElements, 0.);
				std::shared_ptr<MultiComplex> tmp = paras.GetMultiComplex(0);
				for(unsigned int ele=0; ele<nElements; ele++)
					results[ele] = std::norm(tmp->GetValue(ele));

				out = std::shared_ptr<AbsParameter>(new MultiDouble(out->GetName(),results));
			}

			break;
		}//end multi double

		case ParType::DOUBLE:{
			//output double: norm of one double or complex input
			if(!nD && !nC)
				throw std::runtime_error("AbsSquare Strategy: Requested output is DOUBLE but no DOUBLE or COMPLEX was given as input!");
			if(nD){
				std::shared_ptr<DoubleParameter> tmp = paras.GetDoubleParameter(0);
				out = std::shared_ptr<AbsParameter>(new DoubleParameter(out->GetName(),std::norm(tmp->GetValue())));
			}else if(nC){
				std::shared_ptr<ComplexParameter> tmp = paras.GetComplexParameter(0);
				out = std::shared_ptr<AbsParameter>(new DoubleParameter(out->GetName(),std::norm(tmp->GetValue())));
				//          std::cout<<"Strategy AbsSquare (DOUBLE) "<<std::norm(tmp->GetValue())<<std::endl;
			}
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
	Power(ParType in):Strategy(in){
	};
	virtual ~Power(){}

	virtual const std::string to_str() const{
		return "^";
	}

	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out){
		if( checkType != out->type() ) return false;

		unsigned int nMC = paras.GetNMultiComplex();
		unsigned int nMD = paras.GetNMultiDouble();
		unsigned int nC = paras.GetNComplex();
		unsigned int nD = paras.GetNDouble();
		unsigned int nI = paras.GetNInteger();

		if(nMC+nMD+nD+nI==0){
			//TODO: exception no input
			return false;
		}
		//only two double or complex parameter possible
		if( !(nD+nMD+nC+nMC==2) ){
			//TODO: exception wrong input
			return false;
		}

		switch(checkType){

		case ParType::MCOMPLEX:{
			//output multi complex: input must be two multi complex
			if(!(nMC==2)){
				//TODO: exception wrong input
				return false;
			}
			unsigned int nElements = paras.GetMultiComplex(0)->GetNValues();
			//fill MultiDouble parameter
			std::vector<std::complex<double> > results(nElements, std::complex<double>(0.,0.));
			std::shared_ptr<MultiComplex> tmpA = paras.GetMultiComplex(0);
			std::shared_ptr<MultiComplex> tmpB = paras.GetMultiComplex(1);
			for(unsigned int ele=0; ele<nElements; ele++)
				results[ele] = std::pow(tmpA->GetValue(ele),tmpB->GetValue(ele));

			out = std::shared_ptr<AbsParameter>(new MultiComplex(out->GetName(),results));

			break;
		}//end multi complex

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
				results[ele] = std::pow(tmpA->GetValue(ele),tmpB->GetValue(ele));

			out = std::shared_ptr<AbsParameter>(new MultiDouble(out->GetName(),results));

			break;
		}//end multi double

		case ParType::COMPLEX:{
			//output complex: power of two complex input
			if(!(nC==2)){
				//TODO: exception wrong input
				return false;
			}
			std::shared_ptr<ComplexParameter> tmpA = paras.GetComplexParameter(0);
			std::shared_ptr<ComplexParameter> tmpB = paras.GetComplexParameter(1);
			out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(),std::pow(tmpA->GetValue(),tmpB->GetValue())));
			break;
		}//end double

		case ParType::DOUBLE:{
			//output double: power of two double input
			if(!(nD==2)){
				//TODO: exception wrong input
				return false;
			}
			std::shared_ptr<DoubleParameter> tmpA = paras.GetDoubleParameter(0);
			std::shared_ptr<DoubleParameter> tmpB = paras.GetDoubleParameter(1);
			out = std::shared_ptr<AbsParameter>(new DoubleParameter(out->GetName(),std::pow(tmpA->GetValue(),tmpB->GetValue())));
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
