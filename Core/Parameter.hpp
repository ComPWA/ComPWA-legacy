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
//! Implementations for internal parameter.
/*! \class DoubleParameter
 * \class IntegerParameter
 * \class BoolParameter
 * @file Parameter.hpp
 * This class implements some internal container of parameters.
 * A parameter consists of a value with optional bounds and error.
 */

#ifndef _PARAMETER_HPP_
#define _PARAMETER_HPP_

#include <iostream>
#include <string>
#include <sstream>
#include <complex>
#include <stdexcept>
#include <cmath>

#include "Core/AbsParameter.hpp"
#include "Core/Exceptions.hpp"
#include "Core/Logging.hpp"

enum ErrorType {
	SYM = 1, ASYM = 2, LHSCAN = 3, NOTDEF = 0
};

namespace ComPWA {

class MultiComplex: public AbsParameter {

public:

	//! Standard constructor without information
	/*!
	 * Standard constructor with no information provided. Creates parameter
	 * with value 0 but without bounds or an error.
	 * \param inName internal string identifier of this parameter
	 */
	//MultiComplex(std::string inName):AbsParameter(inName, ParType::MDOUBLE){
	//}
	//! Standard constructor with a value
	/*!
	 * Standard constructor with just a value provided. Creates parameter
	 * with given value but without bounds or an error.
	 * \param inName internal string identifier of this parameter
	 * \param values input vector of values of the parameter
	 */
	MultiComplex(std::string inName,
			const std::vector<std::complex<double> >& values) :
				AbsParameter(inName, ParType::MCOMPLEX), val_(values) {
	}

	//! Copy constructor using = operator
	/*!
	 * Simple copy constructor using the = operator. As this operator is not
	 * overloaded in this class, c++ will copy every member variable. As this
	 * is a container class, this should be fine.
	 * \param in input PWAParameter which variables will be copied
	 */
	MultiComplex(const MultiComplex& in) :
		AbsParameter(in.name_, ParType::MCOMPLEX) {
		*this = in;
		//      error_ = std::shared_ptr<ParError<double>>(new ParError<double>(*in.error_));
	}

	//! Empty Destructor
	/*!
	 * There is nothing to destroy :(
	 */
	virtual ~MultiComplex() { /* nothing */
	}

	//! Getter for number of values in this multipar
	virtual const inline double GetNValues() const {
		return val_.size();
	}

	//! Getter for value of parameter
	virtual const inline std::vector<std::complex<double> >& GetValues() const {
		return val_;
	}

	//! Getter for value of parameter
	virtual const inline std::complex<double> GetValue(unsigned int i = 0) const {
		if (i >= val_.size())
			return 0;
		return val_[i];
	}

	//! Getter for FunctionTree support
	virtual const std::complex<double> getNodeValue(unsigned int i = 0) {
		if (i >= val_.size())
			return std::complex<double>();
		return val_[i];
	}

	//! Setter for value of parameter
	virtual void SetValue(const std::complex<double> inVal, unsigned int i = 0) {
		if (i >= val_.size())
			return;
		if (val_[i] == inVal)
			return;
		val_[i] = inVal;
		Notify();
	}

	std::vector<std::complex<double> >::iterator Begin() { return val_.begin(); }
	std::vector<std::complex<double> >::iterator End() { return val_.end(); }

protected:

	virtual std::string TypeName() const {
		return "complex collection";
	}
	//std::string out_; /*!< Output string to print information */
	std::vector<std::complex<double> > val_;/*!< Containers of parameter information */

	//! A protected function which returns an output string for printing
	/*!
	 * This function uses all available information about the parameter
	 * to create a string which will be streamed via the stream operator <<.
	 * \sa operator<<, to_str(), type()
	 */
	virtual std::string make_str() {
		std::stringstream oss;
		oss << name_;
		unsigned int max = val_.size();
		if (max > 5)
			max = 5;    //display only 5 variables
		oss << "\t Val = ";
		for (unsigned int i = 0; i < max - 1; i++)
			oss << val_[i] << ", ";
		oss << val_[max - 1];
		if (max < val_.size())
			oss << " ... ";
		oss << "\t Type = " << TypeName();
		return oss.str();
	}

	//! A protected function which returns an output string for printing
	/*!
	 * This function uses only the value information about the parameter
	 * to create a string which will be streamed via the stream operator <<.
	 * \sa operator<<, make_str()
	 */
	virtual std::string make_val_str() {
		std::stringstream ovs;
		unsigned int max = val_.size();
		if (max > 0) {
			if (max > 3)
				max = 3;    //display only 10 variables
			for (unsigned int i = 0; i < max - 1; i++)
				ovs << val_[i] << ", ";
			ovs << val_[max - 1];
			if (max < val_.size())
				ovs << " ... ";
		}
		return ovs.str();
	}

private:

};

class MultiDouble: public AbsParameter {

public:

	//! Standard constructor without information
	/*!
	 * Standard constructor with no information provided. Creates parameter
	 * with value 0 but without bounds or an error.
	 * \param inName internal string identifier of this parameter
	 */
	//MultiDouble(std::string inName):AbsParameter(inName, ParType::MDOUBLE){
	//}
	//! Standard constructor with a value
	/*!
	 * Standard constructor with just a value provided. Creates parameter
	 * with given value but without bounds or an error.
	 * \param inName internal string identifier of this parameter
	 * \param values input vector of values of the parameter
	 */
	MultiDouble(std::string inName, const std::vector<double>& values) :
		AbsParameter(inName, ParType::MDOUBLE), val_(values) {
	}

	//! Copy constructor using = operator
	/*!
	 * Simple copy constructor using the = operator. As this operator is not
	 * overloaded in this class, c++ will copy every member variable. As this
	 * is a container class, this should be fine.
	 * \param in input PWAParameter which variables will be copied
	 */
	MultiDouble(const MultiDouble& in) :
		AbsParameter(in.name_, ParType::MDOUBLE) {
		*this = in;
		//		error_ = std::shared_ptr<ParError<double>>(new ParError<double>(*in.error_));
	}

	//! Empty Destructor
	/*!
	 * There is nothing to destroy :(
	 */
	virtual ~MultiDouble() { /* nothing */
	}

	//! Getter for number of values in this multipar
	virtual const inline double GetNValues() const {
		return val_.size();
	}

	//! Getter for value of parameter
	virtual const inline std::vector<double>& GetValues() const {
		return val_;
	}

	//! Getter for value of parameter
	virtual const inline double GetValue(unsigned int i = 0) const {
		if (i >= val_.size())
			return 0;
		return val_[i];
	}

	//! Getter for FunctionTree support
	virtual const std::complex<double> getNodeValue(unsigned int i = 0) {
		if (i >= val_.size())
			return std::complex<double>();
		return std::complex<double>(val_[i], 0.);
	}

	//! Setter for value of parameter
	virtual void SetValue(const double inVal, unsigned int i = 0) {
		if (i >= val_.size())
			return;
		if (val_[i] == inVal)
			return;
		val_[i] = inVal;
		Notify();
	}

	std::vector<double>::iterator Begin() { return val_.begin(); }
	std::vector<double>::iterator End() { return val_.end(); }

protected:
	virtual std::string TypeName() const {
		return "double collection";
	}
	//std::string out_; /*!< Output string to print information */
	std::vector<double> val_;/*!< Containers of parameter information */

	//! A protected function which returns an output string for printing
	/*!
	 * This function uses all available information about the parameter
	 * to create a string which will be streamed via the stream operator <<.
	 * \sa operator<<, to_str(), type()
	 */
	virtual std::string make_str() {
		std::stringstream oss;
		oss << name_;
		unsigned int max = val_.size();
		if (max > 5)
			max = 5;    //display only 10 variables
		oss << "\t Val = ";
		for (unsigned int i = 0; i < max - 1; i++)
			oss << val_[i] << ", ";
		oss << val_[max - 1];
		if (max < val_.size())
			oss << " ... ";
		oss << "\t Type = " << TypeName();
		return oss.str();
	}

	//! A protected function which returns an output string for printing
	/*!
	 * This function uses only the value information about the parameter
	 * to create a string which will be streamed via the stream operator <<.
	 * \sa operator<<, make_str()
	 */
	virtual std::string make_val_str() {
		std::stringstream ovs;
		unsigned int max = val_.size();
		if (max > 0) {
			if (max > 5)
				max = 5;    //display only 5 variables
			for (unsigned int i = 0; i < max - 1; i++)
				ovs << val_[i] << ", ";
			ovs << val_[max - 1];
			if (max < val_.size())
				ovs << " ... ";
		}
		return ovs.str();
	}

private:
};

class MultiUnsignedInteger: public AbsParameter {

public:
	//! Standard constructor with a value
	/*!
	 * Standard constructor with just a value provided. Creates parameter
	 * with given value but without bounds or an error.
	 * \param inName internal string identifier of this parameter
	 * \param values input vector of values of the parameter
	 */
	MultiUnsignedInteger(std::string inName,
			const std::vector<unsigned int>& values) :
				AbsParameter(inName, ParType::MUNSIGNEDINTEGER), val_(values) {
	}

	//! Copy constructor using = operator
	/*!
	 * Simple copy constructor using the = operator. As this operator is not
	 * overloaded in this class, c++ will copy every member variable. As this
	 * is a container class, this should be fine.
	 * \param in input PWAParameter which variables will be copied
	 */
	MultiUnsignedInteger(const MultiUnsignedInteger& in) :
		AbsParameter(in.name_, ParType::MUNSIGNEDINTEGER) {
		*this = in;
		//    error_ = std::shared_ptr<ParError<double>>(new ParError<double>(*in.error_));
	}

	//! Empty Destructor
	/*!
	 * There is nothing to destroy :(
	 */
	virtual ~MultiUnsignedInteger() { /* nothing */
	}

	//! Getter for number of values in this multipar
	virtual const inline double GetNValues() const {
		return val_.size();
	}

	//! Getter for value of parameter
	virtual const inline std::vector<unsigned int>& GetValues() const {
		return val_;
	}

	//! Getter for value of parameter
	virtual const inline unsigned int GetValue(unsigned int i = 0) const {
		if (i >= val_.size())
			return 0;
		return val_[i];
	}

	//! Setter for value of parameter
	virtual void SetValue(const unsigned int inVal, unsigned int i = 0) {
		if (i >= val_.size())
			return;
		if (val_[i] == inVal)
			return;
		val_[i] = inVal;
		Notify();
	}

protected:
	virtual std::string TypeName() const {
		return "unsigned int collection";
	}
	//std::string out_; /*!< Output string to print information */
	std::vector<unsigned int> val_;/*!< Containers of parameter information */

	//! A protected function which returns an output string for printing
	/*!
	 * This function uses all available information about the parameter
	 * to create a string which will be streamed via the stream operator <<.
	 * \sa operator<<, to_str(), type()
	 */
	virtual std::string make_str() {
		std::stringstream oss;
		oss << name_;
		unsigned int max = val_.size();
		if (max > 5)
			max = 5;    //display only 10 variables
		oss << "\t Val = ";
		for (unsigned int i = 0; i < max - 1; i++)
			oss << val_[i] << ", ";
		oss << val_[max - 1];
		if (max < val_.size())
			oss << " ... ";
		oss << "\t Type = " << TypeName();
		return oss.str();
	}

	//! A protected function which returns an output string for printing
	/*!
	 * This function uses only the value information about the parameter
	 * to create a string which will be streamed via the stream operator <<.
	 * \sa operator<<, make_str()
	 */
	virtual std::string make_val_str() {
		std::stringstream ovs;
		unsigned int max = val_.size();
		if (max > 5)
			max = 5;    //display only 5 variables
		for (unsigned int i = 0; i < max - 1; i++)
			ovs << val_[i] << ", ";
		ovs << val_[max - 1];
		if (max < val_.size())
			ovs << " ... ";
		return ovs.str();
	}

private:
};

class ComplexParameter: public AbsParameter {

public:

	//! Standard constructor without information
	/*!
	 * Standard constructor with just a name provided. Creates parameter
	 * with value 0 but without bounds or an error.
	 * \param inName internal string identifier of this parameter
	 */
	ComplexParameter(std::string inName) :
		AbsParameter(inName, ParType::COMPLEX), val_(0., 0.), min_(0., 0.), max_(
				0., 0.), err_(0., 0.) {
		bounds_ = usebounds_ = hasError_ = fixed_ = false;
	}

	//! Standard constructor with a value
	/*!
	 * Standard constructor with just value and name provided. Creates parameter
	 * with given value but without bounds or an error.
	 * \param inName internal string identifier of this parameter
	 * \param value input value of the parameter
	 */
	ComplexParameter(std::string inName, const std::complex<double> value) :
		AbsParameter(inName, ParType::COMPLEX), val_(value), min_(0, 0), max_(0,
				0), err_(0, 0) {
		bounds_ = usebounds_ = hasError_ = fixed_ = false;
	}

	//! Standard constructor with value and error
	/*!
	 * Standard constructor with value and error provided. Creates parameter
	 * with given value and error but without bounds.
	 * \param inName internal string identifier of this parameter
	 * \param value input value of the parameter
	 * \param error input error of the parameter
	 */
	ComplexParameter(std::string inName, const std::complex<double> value,
			const std::complex<double> error) :
				AbsParameter(inName, ParType::COMPLEX), val_(value), min_(0, 0), max_(0,
						0), err_(error) {
		bounds_ = usebounds_ = fixed_ = false;
		hasError_ = true;
	}

	//! Standard constructor with value and bounds
	/*!
	 * Standard constructor with value and bounds provided. Creates parameter
	 * with given value and bounds but without error. If a check for valid
	 * bounds fails, just the value is used.
	 * \param inName internal string identifier of this parameter
	 * \param value input value of the parameter
	 * \param min input lower bound
	 * \param max input upper bound
	 * \sa check_bounds()
	 */
	ComplexParameter(std::string inName, const std::complex<double> value,
			const std::complex<double> min, const std::complex<double> max) :
				AbsParameter(inName, ParType::COMPLEX), val_(value), min_(0, 0), max_(0,
						0), err_(0, 0) {
		bounds_ = usebounds_ = hasError_ = fixed_ = false;
		if (check_bounds(min, max)) {
			min_ = min;
			max_ = max;
			bounds_ = true;
		}
	}

	//! Standard constructor with value, bounds and error
	/*!
	 * Standard constructor with value, bounds and error provided. Creates
	 * parameter with the given information. If a check for valid bounds
	 * fails, just value and error are used.
	 * \param inName internal string identifier of this parameter
	 * \param value input value of the parameter
	 * \param min input lower bound
	 * \param max input upper bound
	 * \param error input error of the parameter
	 * \sa check_bounds()
	 */
	ComplexParameter(std::string inName, const std::complex<double> value,
			const std::complex<double> min, const std::complex<double> max,
			const std::complex<double> error) :
				AbsParameter(inName, ParType::COMPLEX), val_(value), min_(0, 0), max_(0,
						0), err_(error) {
		bounds_ = usebounds_ = fixed_ = false;
		hasError_ = true;
		if (check_bounds(min, max)) {
			min_ = min;
			max_ = max;
			bounds_ = true;
		}
	}

	//! Copy constructor using = operator
	/*!
	 * Simple copy constructor using the = operator. As this operator is not
	 * overloaded in this class, c++ will copy every member variable. As this
	 * is a container class, this should be fine.
	 * \param in input PWAParameter which variables will be copied
	 */
	ComplexParameter(const ComplexParameter& in) :
		AbsParameter(in.name_, ParType::COMPLEX) {
		*this = in;
	}

	//! Empty Destructor
	/*!
	 * There is nothing to destroy :(
	 */
	virtual ~ComplexParameter() { /* nothing */
	}

	//! Check if parameter has bounds
	virtual const inline bool HasBounds() const {
		return bounds_;
	}
	//! Check if bounds should be used
	virtual const inline bool UseBounds() const {
		if (bounds_)
			return usebounds_;
		return false;
	}
	//! Check if parameter has an error
	virtual const inline bool HasError() const {
		return hasError_;
	}
	//! Check if parameter is fixed
	virtual const inline bool IsFixed() const {
		return fixed_;
	}

	//! Getter for value of parameter
	virtual const inline std::complex<double> GetValue() const {
		return val_;
	}
	//! Getter for lower bound of parameter
	virtual const inline std::complex<double> GetMinValue() const {
		return min_;
	}
	//! Getter for upper bound of parameter
	virtual const inline std::complex<double> GetMaxValue() const {
		return max_;
	}
	//! Getter for error of parameter
	virtual const inline std::complex<double> GetError() const {
		return err_;
	}

	//! Getter for FunctionTree support
	//virtual const std::complex<double> getNodeValue(){
	//  return val_;
	//}

	//! Setter for value of parameter
	virtual void SetValue(const std::complex<double> inVal) {
		if (fixed_) {
			throw ParameterFixed();
			return;
		}
		if (val_ == inVal)
			return;
		val_ = inVal;
		Notify();
	}
	//! Setter for error of parameter
	virtual void SetError(const std::complex<double> inErr) {
		err_ = inErr;
		hasError_ = true;
	}    //make_str();}
	//! Setter for bounds of parameter
	virtual const bool SetMinMax(const std::complex<double> inMin,
			const std::complex<double> inMax) {
		bool valid = check_bounds(inMin, inMax);
		if (valid) {
			min_ = inMin;
			max_ = inMax;
			bounds_ = true;
		}
		return valid;
	}

	//! Setter for lower bound
	/*!
	 * Setter for lower bound of the parameter. If a check for valid bounds
	 * fails, it returns false and nothing changes. This means if the lower
	 * bound is invalid the parameter maintains its old bounds if it had some.
	 * \param min input lower bound
	 * \return bool if successful (re)set lower bound
	 * \sa check_bounds()
	 */
	virtual const bool SetMinValue(const std::complex<double> min) {
		bool valid = check_bounds(min, max_);
		if (valid) {
			min_ = min;
			bounds_ = true;
		}
		return valid;
	}

	//! Setter for upper bound
	/*!
	 * Setter for upper bound of the parameter. If a check for valid bounds
	 * fails, it returns false and nothing changes. This means if the upper
	 * bound is invalid the parameter maintains its old bounds if it had some.
	 * \param max input upper bound
	 * \return bool if successful (re)set upper bound
	 * \sa check_bounds()
	 */
	virtual const bool SetMaxValue(const std::complex<double> max) {
		bool valid = check_bounds(min_, max);
		if (valid) {
			max_ = max;
			bounds_ = true;
		}
		return valid;
	}

	//! Set if bounds should be used
	virtual const inline void UseBounds(const bool use) {
		usebounds_ = use;
	}
	//! Call to fix parameter
	virtual const inline void SetParameterFixed() {
		fixed_ = true;
	}
	//! Call to free parameter
	virtual const inline void SetParameterFree() {
		fixed_ = false;
	}
	//! Set parameter free or fixed
	virtual const inline void FixParameter(const bool fixed) {
		fixed_ = fixed;
	}

	//! A public function returning a string naming its type
	/*!
	 * This function is used to get the type of the implementation of this
	 * general parameter interface.
	 * \sa operator<<, to_str(), make_str()
	 */
	virtual std::string TypeName() const {
		return "complex double";
	}

protected:
	bool bounds_; /*!< Are valid bounds defined for this parameter? */
	bool hasError_; /*!< Is an error defined for this parameter? */
	bool usebounds_; /*!< Do you want to restrict your parameter? */
	bool fixed_; /*!< Do you want to keep parameter fixed? */
	std::complex<double> val_, min_, max_, err_; /*!< Containers of parameter information */

	//! A protected function to check if bounds are valid
	/*!
	 * This function checks if the bounds of the parameter are valid:
	 * Upper bound should be larger then lower bound and the value
	 * should be inside of the bounds.
	 * \param max upper bound to check
	 * \param min lower bound to check
	 * \return bool if bounds are valid
	 * \sa Parameter(const double value, const double min, const double max)
	 * \sa Parameter(const double value, const double min, const double max, const double error)
	 * \sa SetMinMax(), SetMinValue(), SetMaxValue()
	 */
	bool check_bounds(const std::complex<double> min,
			const std::complex<double> max) {
		if ((max.real() > min.real()) && (max.real() >= val_.real())
				&& (min.real() <= val_.real()) && (max.imag() > min.imag())
				&& (max.imag() >= val_.imag()) && (min.imag() <= val_.imag()))
			return true;
		return false;
	}

	//! A protected function which returns an output string for printing
	/*!
	 * This function uses all available information about the parameter
	 * to create a string which will be streamed via the stream operator <<.
	 * \sa operator<<, to_str(), type()
	 */
	virtual std::string make_str() {
		std::stringstream oss;
		oss << name_;
		oss << "\t Val = " << val_;
		if (bounds_)
			oss << "\t  Min-Max = " << min_ << " to " << max_;
		if (hasError_)
			oss << "\t  Err = " << err_;
		oss << "\t Type = " << TypeName();
		return oss.str();
	}

	//! A protected function which returns an output string for printing
	/*!
	 * This function uses only the value information about the parameter
	 * to create a string which will be streamed via the stream operator <<.
	 * \sa operator<<, make_str()
	 */
	virtual std::string make_val_str() {
		std::stringstream ovs;
		ovs << val_;
		return ovs.str();
	}

};

class DoubleParameter: public AbsParameter {

public:
	//! Standard constructor without information
	/*!
	 * Standard constructor with no information provided. Creates parameter
	 * with value 0 but without bounds or an error.
	 * \param inName internal string identifier of this parameter
	 */
	DoubleParameter(std::string inName="") :
		AbsParameter(inName, ParType::DOUBLE), fixed_(0),val_(0),min_(0), max_(0)
{
		SetError(0);
		bounds_= usebounds_ = false;
}

	//! Standard constructor with a value
	/*!
	 * Standard constructor with just a value provided. Creates parameter
	 * with given value but without bounds or an error.
	 * \param inName internal string identifier of this parameter
	 * \param value input value of the parameter
	 */
	DoubleParameter(std::string inName, const double value) :
		AbsParameter(inName, ParType::DOUBLE), fixed_(0),val_(value), min_(0),max_(0)
	{
		SetError(0);
		bounds_= usebounds_ = false;
	}

	//! Standard constructor with value and error
	/*!
	 * Standard constructor with value and error provided. Creates parameter
	 * with given value and error but without bounds.
	 * \param inName internal string identifier of this parameter
	 * \param value input value of the parameter
	 * \param error input error of the parameter
	 */
	DoubleParameter(std::string inName, const double value, const double error)
	:AbsParameter(inName, ParType::DOUBLE),fixed_(0),val_(value),min_(0),max_(0),
	 errorType(ErrorType::NOTDEF)
	{
		SetError(error);
		bounds_= usebounds_ = false;
	}

	//! Standard constructor with value and bounds
	/*!
	 * Standard constructor with value and bounds provided. Creates parameter
	 * with given value and bounds but without error. If a check for valid
	 * bounds fails, just the value is used.
	 * \param inName internal string identifier of this parameter
	 * \param value input value of the parameter
	 * \param min input lower bound
	 * \param max input upper bound
	 * \sa check_bounds()
	 */
	DoubleParameter(std::string inName, const double value,
			const double min, const double max) :
				AbsParameter(inName, ParType::DOUBLE),
				fixed_(0),val_(value),min_(0),max_(0)
	{
		SetError(0);
		bounds_= usebounds_ = false;
		SetMinMax(min,max);
	}

	//! Standard constructor with value, bounds and error
	/*!
	 * Standard constructor with value, bounds and error provided. Creates
	 * parameter with the given information. If a check for valid bounds
	 * fails, just value and error are used.
	 * \param inName internal string identifier of this parameter
	 * \param value input value of the parameter
	 * \param min input lower bound
	 * \param max input upper bound
	 * \param error input error of the parameter
	 * \sa check_bounds()
	 */
	DoubleParameter(std::string inName, const double value,
			const double min, const double max, const double error) :
				AbsParameter(inName, ParType::DOUBLE),
				fixed_(0),val_(value),min_(0),max_(0)
	{
		SetError(error);
		bounds_= usebounds_ = false;
		SetMinMax(min,max);
	}
	DoubleParameter(const DoubleParameter& in) :
		AbsParameter(in.name_, ParType::DOUBLE)
	{
		*this = in;
	}
	//! Empty Destructor
	virtual ~DoubleParameter() { /* nothing */	}

	//! Operator for conversion to double
	operator double() const{ return val_;};

	//! Check if parameter has bounds
	virtual inline bool HasBounds() const {return bounds_;}
	//! Check if bounds should be used
	virtual inline bool UseBounds() const {
		if(bounds_) return usebounds_;
		return false;
	}
	//! Check if parameter is fixed
	virtual inline bool IsFixed() const {return fixed_;}
	//! Set if bounds should be used
	virtual inline void SetUseBounds(const bool use) {usebounds_=use;}
	//! Call to fix parameter
	virtual inline void SetParameterFixed() {fixed_=true;}
	//! Call to free parameter
	virtual inline void SetParameterFree() {fixed_=false;}
	//! Set parameter free or fixed
	virtual inline void FixParameter(const bool fixed) {fixed_=fixed;}
	/*! Update member variables from other DoubleParameter
	 * Do to the Observer pattern we can't use a copy constructor.
	 * Therefore we use this workaround. The function ignores if parameter
	 * is fixed!
	 */
	virtual void UpdateParameter( std::shared_ptr<DoubleParameter> newPar ){
		//copy bounds
		if(newPar->HasBounds()){
			try{
				SetMinMax(newPar->GetMinValue(), newPar->GetMaxValue());
			} catch (std::exception& ex){
				//ignore if bound can not be set
			}
			SetUseBounds(newPar->UseBounds());
		} else
			bounds_ = usebounds_ = 0;

		bool isFix = newPar->IsFixed();
		//copy value
		FixParameter(0); //we ignore here if parameter is fixed
		SetValue(newPar->GetValue());

		//Check if bounds are valid and value is within bounds
		if( !check_bounds(GetMinValue(),GetMaxValue()) )
			throw std::runtime_error("DoubleParameter::UpdateParameter() | "
					"Bounds not valid for parameter "+GetName()+": "
					+std::to_string(GetValue())
		+" ["+std::to_string((long double)GetMinValue())+";"
		+std::to_string((long double)GetMaxValue())+"]!");

		//copy error
		if(newPar->GetErrorType()==ErrorType::SYM)
			SetError(newPar->GetError());
		else if(newPar->GetErrorType()==ErrorType::ASYM)
			SetError(newPar->GetErrorLow(),newPar->GetErrorHigh());
		else
			SetErrorType(ErrorType::NOTDEF);

		//copy fix parameter
		FixParameter(isFix);
		return;
	}

	//====== PARAMETER VALUE ========
	//! Getter for value of parameter
	virtual inline double GetValue() const {return val_;}
	//! Getter for value of parameter
	virtual inline double GetRoundedValue() const {
		return val_;
	}
	//! Getter for lower bound of parameter
	virtual inline double GetMinValue() const {return min_;}
	//! Getter for upper bound of parameter
	virtual inline double GetMaxValue() const {return max_;}
	//! Getter for FunctionTree support
	virtual std::complex<double> getNodeValue(){
		return std::complex<double>(val_,0.);
	}
	//! Setter for value of parameter
	virtual void SetValue(const double inVal) {
		if(fixed_)
			throw ParameterFixed("DoubleParameter::SetValue() | Parameter "+GetName()+" is fixed!");
		/*Call notify only if value has changed! Otherwise tree is
		 * recalculated also in case where current parameter is not changed
		 */
		//if(std::fabs(val_-inVal) < 0.0000001) return;
		if(val_==inVal) return;

		if(usebounds_ && (inVal < GetMinValue() || inVal > GetMaxValue()) )
			throw ParameterOutOfBound("DoubleParameter::SetValue() | "
					"Parameter "+GetName()+" not within bounds: val="
					+std::to_string(inVal)
		+" ["+std::to_string((long double)GetMinValue())+";"
		+std::to_string((long double)GetMaxValue())+"]!");

		val_ = inVal;
		Notify();
	}

	//! Setter for bounds of parameter
	virtual void SetMinMax(const double min, const double max){
		try{
			SetMinValue(min);
		} catch (ParameterOutOfBound& ex) { }
		try{
			SetMaxValue(max);
		} catch (ParameterOutOfBound& ex) { }
		if(!check_bounds(min_,max_))
			throw std::runtime_error("DoubleParameter::SetMinMaxValue() | "
					"Bounds not valid for parameter "+GetName()+": "
					+std::to_string(GetValue())
		+" ["+std::to_string((long double)min_)+";"
		+std::to_string((long double)max_)+"]!");
		bounds_ = usebounds_ = true;
	}
	/*! Setter for lower bound
	 * Setter for lower bound of the parameter. If a check for valid bounds
	 * fails, it returns false and nothing changes. This means if the lower
	 * bound is invalid the parameter maintains its old bounds if it had some.
	 * \param min input lower bound
	 * \return bool if successful (re)set lower bound
	 * \sa check_bounds()
	 */
	virtual void SetMinValue(const double min) {
		min_ = min;
		if( !check_bounds(min_, max_) )
			throw ParameterOutOfBound("DoubleParameter::SetMinValue() | "
					"Boundary not valid: ["
					+std::to_string(GetMinValue())+", "
					+std::to_string(GetMaxValue())+"]!"
			);
		bounds_ = usebounds_ = true;
	}
	/*! Setter for upper bound
	 * Setter for upper bound of the parameter. If a check for valid bounds
	 * fails, it returns false and nothing changes. This means if the upper
	 * bound is invalid the parameter maintains its old bounds if it had some.
	 * \param max input upper bound
	 * \return bool if successful (re)set upper bound
	 * \sa check_bounds()
	 */
	virtual void SetMaxValue(const double max) {
		max_ = max;
		if( !check_bounds(min_, max_) )
			throw ParameterOutOfBound("DoubleParameter::SetMaxValue() | "
					"Boundary not valid: ["
					+std::to_string(GetMinValue())+", "
					+std::to_string(GetMaxValue())+"]!"
			);
		bounds_ = usebounds_ = true;
	}
	//====== PARAMETER ERROR ========
	//! Check if parameter has an error
	virtual inline bool HasError() const {
		if(GetErrorType()==ErrorType::NOTDEF) return 0;
		else return 1;
	}
	//! Getter for type of parameter error
	virtual ErrorType GetErrorType() const {return errorType;}
	//! Getter for parameter error. In case of asymmetric errors the average error is returned.
	virtual double GetError() const {
		if(!HasError())
			throw std::runtime_error("DoubleParameter::GetError() | "
					"Parameter "+name_+" has no errors defined!");
		return (GetErrorHigh()+GetErrorLow())/2;
	}
	//! Get rounded parameter error. In case of asymmetric errors the average error is returned.
	virtual double GetRoundedError() const {
		return GetError();
	}
	//! Getter for upper error of parameter
	virtual double GetErrorHigh() const {
		if(!HasError())
			throw std::runtime_error("DoubleParameter::GetError() | "
					"Parameter "+name_+" has no errors defined!");
		//		if(GetErrorType()==ErrorType::SYM){
		//			BOOST_LOG_TRIVIAL(info) << "DoubleParameter::GetErrorHigh() | Parameter "<<name_
		//					<<" has no asymmetric errors! Returning symmetric error";
		//			return GetError();
		//		}
		return errorHigh;
	}
	//! Getter for lower error of parameter
	virtual double GetErrorLow() const {
		GetName();
		if(!HasError())
			throw std::runtime_error("DoubleParameter::GetError() | "
					"Parameter "+name_+" has no errors defined!");
		//		if(GetErrorType()==ErrorType::SYM){
		//			BOOST_LOG_TRIVIAL(info) << "DoubleParameter::GetErrorHigh() | Parameter "<<name_
		//					<<" has no assymetric errors! Returning symmetric error";
		//			return GetError();
		//		}
		if(!HasError())
			throw std::runtime_error("DoubleParameter::GetError() | "
					"Parameter "+name_+" has no errors defined!");
		return errorLow;
	}
	//! Setter for low/high error of parameter
	virtual void SetError(double errLow, double errHigh) {
		//if(fixed_)
		//	throw ParameterFixed("DoubleParameter::SetError(double, double) | Parameter "+GetName()+" is fixed!");
		SetErrorType(ErrorType::ASYM);
		SetErrorHigh(errHigh);
		SetErrorLow(errLow);
	}
	//! Setter for error of parameter
	virtual void SetError(double err) {
		//if(fixed_)
		//	throw ParameterFixed("DoubleParameter::SetError(double) | Parameter "+GetName()+" is fixed!");
		SetErrorType(ErrorType::SYM);
		SetErrorHigh(err);
		SetErrorLow(err);
	}

protected:
	//	std::string out_; /*!< Output string to print information */
	/*! A public function returning a string naming its type
	 * This function is used to get the type of the implementation of this
	 * general parameter interface.
	 * \sa operator<<, to_str(), make_str()
	 */
	virtual std::string TypeName() const { return "double"; }
	bool bounds_; /*!< Are valid bounds defined for this parameter? */
	bool usebounds_; /*!< Do you want to restrict your parameter? */
	bool fixed_; /*!< Do you want to keep parameter fixed? */
	double val_, min_, max_;/*!< Current value, bounds*/
	//! error type
	ErrorType errorType;
	//! lower parameter error
	double errorLow;
	//! upper parameter error
	double errorHigh;
	//! Setter for upper error of parameter
	virtual void SetErrorHigh(double errHigh) { errorHigh=errHigh;}
	//! Setter for lower error of parameter
	virtual void SetErrorLow(double errLow) { errorLow=std::fabs(errLow);}
	//! Setter for type of parameter error
	virtual void SetErrorType(ErrorType t) { errorType=t; }

	//! A protected function to check if bounds are valid
	/*!
	 * This function checks if the bounds of the parameter are valid:
	 * Upper bound should be larger then lower bound and the value
	 * should be inside of the bounds.
	 * \param max upper bound to check
	 * \param min lower bound to check
	 * \return bool if bounds are valid
	 * \sa Parameter(const double value, const double min, const double max)
	 * \sa Parameter(const double value, const double min, const double max, const double error)
	 * \sa SetMinMax(), SetMinValue(), SetMaxValue()
	 */
	bool check_bounds(const double min, const double max){
		if( (max > min) && (max >= val_) && (min <= val_))
			return true;
		return false;
	}

	//! A protected function which returns an output string for printing
	/*!
	 * This function uses all available information about the parameter
	 * to create a string which will be streamed via the stream operator <<.
	 * \sa operator<<, to_str(), type()
	 */
	virtual std::string make_str() {
		std::stringstream oss;
		oss << name_;
		oss << "\t Val = " << val_;
		if(errorLow!=errorHigh)
			oss << " (+"<< errorHigh<<" -"<<errorLow<<")";
		else if ( errorLow!=0 )
			oss << " (+-"<< errorLow<<")";
		if(bounds_)
			oss << "\t  [" << min_ << " ; " << max_<<"]";
		oss << " fix? "<<IsFixed();
		oss << "\t " << TypeName();
		return oss.str();
	}

	//! A protected function which returns an output string for printing
	/*!
	 * This function uses only the value information about the parameter
	 * to create a string which will be streamed via the stream operator <<.
	 * \sa operator<<, make_str()
	 */
	virtual std::string make_val_str() {
		std::stringstream ovs;
		ovs << val_;
		return ovs.str();
	}

private:
	friend class boost::serialization::access;
	template<class archive>
	void serialize(archive& ar, const unsigned int version)
	{
		using namespace boost::serialization;
		ar & boost::serialization::make_nvp(
				"AbsParameter",
				boost::serialization::base_object<AbsParameter>(*this)
		);
		//		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbsParameter);  //serialize base class
		ar & make_nvp("bounds",bounds_);
		ar & make_nvp("usebounds",usebounds_);
		ar & make_nvp("isFixed",fixed_);
		ar & make_nvp("value",val_);
		ar & make_nvp("min_value",min_);
		ar & make_nvp("max_value",max_);
		try{
			ar & make_nvp("errorType",errorType);
			ar & make_nvp("errorLow",errorLow);
			ar & make_nvp("errorHigh",errorHigh);
		} catch (...) {
			errorLow = 0;
			errorHigh = 0;
			errorType = ErrorType::SYM;
		}
	}
};

class IntegerParameter: public AbsParameter {

public:

	//! Standard constructor without information
	/*!
	 * Standard constructor with no information provided. Creates parameter
	 * with value 0 but without bounds or an error.
	 * \param inName internal string identifier of this parameter
	 */
	IntegerParameter(std::string inName) :
		AbsParameter(inName, ParType::INTEGER), val_(0), min_(0), max_(0), err_(0) {
		bounds_ = usebounds_ = hasError_ = fixed_ = false;
	}

	//! Standard constructor with a value
	/*!
	 * Standard constructor with just a value provided. Creates parameter
	 * with given value but without bounds or an error.
	 * \param inName internal string identifier of this parameter
	 * \param value input value of the parameter
	 */
	IntegerParameter(std::string inName, const int value) :
		AbsParameter(inName, ParType::INTEGER), val_(value), min_(0), max_(0), err_(
				0) {
		bounds_ = usebounds_ = hasError_ = fixed_ = false;
	}

	//! Standard constructor with value and error
	/*!
	 * Standard constructor with value and error provided. Creates parameter
	 * with given value and error but without bounds.
	 * \param inName internal string identifier of this parameter
	 * \param value input value of the parameter
	 * \param error input error of the parameter
	 */
	IntegerParameter(std::string inName, const int value, const int error) :
		AbsParameter(inName, ParType::INTEGER), val_(value), min_(0), max_(0), err_(
				error) {
		bounds_ = usebounds_ = fixed_ = false;
		hasError_ = true;
	}

	//! Standard constructor with value and bounds
	/*!
	 * Standard constructor with value and bounds provided. Creates parameter
	 * with given value and bounds but without error. If a check for valid
	 * bounds fails, just the value is used.
	 * \param inName internal string identifier of this parameter
	 * \param value input value of the parameter
	 * \param min input lower bound
	 * \param max input upper bound
	 * \sa check_bounds()
	 */
	IntegerParameter(std::string inName, const int value, const int min,
			const int max) :
				AbsParameter(inName, ParType::INTEGER), val_(value), min_(0), max_(0), err_(
						0) {
		bounds_ = usebounds_ = hasError_ = fixed_ = false;
		if (check_bounds(min, max)) {
			min_ = min;
			max_ = max;
			bounds_ = true;
		}
	}

	//! Standard constructor with value, bounds and error
	/*!
	 * Standard constructor with value, bounds and error provided. Creates
	 * parameter with the given information. If a check for valid bounds
	 * fails, just value and error are used.
	 * \param inName internal string identifier of this parameter
	 * \param value input value of the parameter
	 * \param min input lower bound
	 * \param max input upper bound
	 * \param error input error of the parameter
	 * \sa check_bounds()
	 */
	IntegerParameter(std::string inName, const int value, const int min,
			const int max, const int error) :
				AbsParameter(inName, ParType::INTEGER), val_(value), min_(0), max_(0), err_(
						error) {
		bounds_ = usebounds_ = fixed_ = false;
		hasError_ = true;
		if (check_bounds(min, max)) {
			min_ = min;
			max_ = max;
			bounds_ = true;
		}
	}

	//! Copy constructor using = operator
	/*!
	 * Simple copy constructor using the = operator. As this operator is not
	 * overloaded in this class, c++ will copy every member variable. As this
	 * is a container class, this should be fine.
	 * \param in input PWAParameter which variables will be copied
	 */
	IntegerParameter(const IntegerParameter& in) :
		AbsParameter(in.name_, ParType::INTEGER) {
		*this = in;
	}

	//! Empty Destructor
	/*!
	 * There is nothing to destroy :(
	 */
	virtual ~IntegerParameter() { /* nothing */
	}

	//! Check if parameter has bounds
	virtual const inline bool HasBounds() const {
		return bounds_;
	}
	//! Check if bounds should be used
	virtual const inline bool UseBounds() const {
		if (bounds_)
			return usebounds_;
		return false;
	}
	//! Check if parameter has an error
	virtual const inline bool HasError() const {
		return hasError_;
	}
	//! Check if parameter is fixed
	virtual const inline bool IsFixed() const {
		return fixed_;
	}

	//! Getter for value of parameter
	virtual const inline int GetValue() const {
		return val_;
	}
	//! Getter for lower bound of parameter
	virtual const inline int GetMinValue() const {
		return min_;
	}
	//! Getter for upper bound of parameter
	virtual const inline int GetMaxValue() const {
		return max_;
	}
	//! Getter for error of parameter
	virtual const inline int GetError() const {
		return err_;
	}

	//! Getter for FunctionTree support
	virtual const std::complex<double> getNodeValue() {
		return std::complex<double>((double) val_, 0.);
	}

	//! Setter for value of parameter
	virtual void SetValue(const int inVal) {
		if (fixed_) {
			throw ParameterFixed();
			return;
		}
		if (val_ == inVal)
			return;
		val_ = inVal;
		Notify();
	}
	//! Setter for error of parameter
	virtual void SetError(const int inErr) {
		err_ = inErr;
		hasError_ = true;
	}    //make_str();}
	//! Setter for bounds of parameter
	virtual bool SetMinMax(const int inMin, const int inMax) {
		bool valid = check_bounds(inMin, inMax);
		if (valid) {
			min_ = inMin;
			max_ = inMax;
			bounds_ = true;
		}
		return valid;
	}

	//! Setter for lower bound
	/*!
	 * Setter for lower bound of the parameter. If a check for valid bounds
	 * fails, it returns false and nothing changes. This means if the lower
	 * bound is invalid the parameter maintains its old bounds if it had some.
	 * \param min input lower bound
	 * \return bool if successful (re)set lower bound
	 * \sa check_bounds()
	 */
	virtual bool SetMinValue(const int min) {
		bool valid = check_bounds(min, max_);
		if (valid) {
			min_ = min;
			bounds_ = true;
		}
		return valid;
	}

	//! Setter for upper bound
	/*!
	 * Setter for upper bound of the parameter. If a check for valid bounds
	 * fails, it returns false and nothing changes. This means if the upper
	 * bound is invalid the parameter maintains its old bounds if it had some.
	 * \param max input upper bound
	 * \return bool if successful (re)set upper bound
	 * \sa check_bounds()
	 */
	virtual bool SetMaxValue(const int max) {
		bool valid = check_bounds(min_, max);
		if (valid) {
			max_ = max;
			bounds_ = true;
		}
		return valid;
	}

	//! Set if bounds should be used
	virtual inline void UseBounds(const bool use) {
		usebounds_ = use;
	}
	//! Call to fix parameter
	virtual inline void SetParameterFixed() {
		fixed_ = true;
	}
	//! Call to free parameter
	virtual inline void SetParameterFree() {
		fixed_ = false;
	}
	//! Set parameter free or fixed
	virtual inline void FixParameter(const bool fixed) {
		fixed_ = fixed;
	}

	//! A public function returning a string naming its type
	/*!
	 * This function is used to get the type of the implementation of this
	 * general parameter interface.
	 * \sa operator<<, to_str(), make_str()
	 */
	operator int() const {
		return val_;
	}
	;

	virtual std::string TypeName() const {
		return "integer";
	}

protected:
	//std::string out_; /*!< Output string to print information */
	bool bounds_; /*!< Are valid bounds defined for this parameter? */
	bool hasError_; /*!< Is an error defined for this parameter? */
	bool usebounds_; /*!< Do you want to restrict your parameter? */
	bool fixed_; /*!< Do you want to keep parameter fixed? */
	int val_, min_, max_, err_; /*!< Containers of parameter information */

	//! A protected function to check if bounds are valid
	/*!
	 * This function checks if the bounds of the parameter are valid:
	 * Upper bound should be larger then lower bound and the value
	 * should be inside of the bounds.
	 * \param max upper bound to check
	 * \param min lower bound to check
	 * \return bool if bounds are valid
	 * \sa Parameter(const double value, const double min, const double max)
	 * \sa Parameter(const double value, const double min, const double max, const double error)
	 * \sa SetMinMax(), SetMinValue(), SetMaxValue()
	 */
	bool check_bounds(const int min, const int max) {
		if ((max > min) && (max >= val_) && (min <= val_))
			return true;
		return false;
	}

	//! A protected function which returns an output string for printing
	/*!
	 * This function uses all available information about the parameter
	 * to create a string which will be streamed via the stream operator <<.
	 * \sa operator<<, to_str(), type()
	 */
	virtual std::string make_str() {
		std::stringstream oss;
		oss << name_;
		oss << "\t Val = " << val_;
		if (bounds_)
			oss << "\t  Min-Max = " << min_ << " to " << max_;
		if (hasError_)
			oss << "\t  Err = " << err_;
		oss << "\t Type = " << TypeName();
		return oss.str();
	}

	//! A protected function which returns an output string for printing
	/*!
	 * This function uses only the value information about the parameter
	 * to create a string which will be streamed via the stream operator <<.
	 * \sa operator<<, make_str()
	 */
	virtual std::string make_val_str() {
		std::stringstream ovs;
		ovs << val_;
		return ovs.str();
	}

};

class BoolParameter: public AbsParameter {

public:

	//! Standard constructor without information
	/*!
	 * Standard constructor with no information provided. Creates parameter
	 * with value 0 but without bounds or an error.
	 * \param inName internal string identifier of this parameter
	 */
	BoolParameter(std::string inName) :
		AbsParameter(inName, ParType::BOOL), val_(0), err_(0) {
		hasError_ = fixed_ = usebounds_ = false;
	}

	//! Standard constructor with a value
	/*!
	 * Standard constructor with just a value provided. Creates parameter
	 * with given value but without bounds or an error.
	 * \param inName internal string identifier of this parameter
	 * \param value input value of the parameter
	 */
	BoolParameter(std::string inName, const bool value) :
		AbsParameter(inName, ParType::BOOL), val_(value), err_(0) {
		hasError_ = fixed_ = usebounds_ = false;
	}

	//! Standard constructor with value and error
	/*!
	 * Standard constructor with value and error provided. Creates parameter
	 * with given value and error but without bounds.
	 * \param inName internal string identifier of this parameter
	 * \param value input value of the parameter
	 * \param error input error of the parameter
	 */
	BoolParameter(std::string inName, const bool value, const bool error) :
		AbsParameter(inName, ParType::BOOL), val_(value), err_(error) {
		usebounds_ = false;
		fixed_ = false;
		hasError_ = true;
	}

	//! Copy constructor using = operator
	/*!
	 * Simple copy constructor using the = operator. As this operator is not
	 * overloaded in this class, c++ will copy every member variable. As this
	 * is a container class, this should be fine.
	 * \param in input PWAParameter which variables will be copied
	 */
	BoolParameter(const BoolParameter& in) :
		AbsParameter(in.name_, ParType::BOOL) {
		*this = in;
	}

	//! Empty Destructor
	/*!
	 * There is nothing to destroy :(
	 */
	virtual ~BoolParameter() { /* nothing */
	}

	//! Check if parameter has an error
	virtual const inline bool HasError() const {
		return hasError_;
	}
	//! Check if parameter is fixed
	virtual const inline bool IsFixed() const {
		return fixed_;
	}

	//! Getter for value of parameter
	virtual const inline bool GetValue() const {
		return val_;
	}
	//! Getter for error of parameter
	virtual const inline bool GetError() const {
		return err_;
	}

	//! Getter for FunctionTree support
	virtual const std::complex<double> getNodeValue() {
		return std::complex<double>((double) val_, 0.);
	}

	//! Setter for value of parameter
	virtual void SetValue(const bool inVal) {
		if (fixed_) {
			throw ParameterFixed();
			return;
		}
		if (val_ == inVal)
			return;
		val_ = inVal;
		Notify();
	}
	//! Setter for error of parameter
	virtual void SetError(const bool inErr) {
		err_ = inErr;
		hasError_ = true;
	}    //make_str();}

	//! Call to fix parameter
	virtual const inline void SetParameterFixed() {
		fixed_ = true;
	}
	//! Call to free parameter
	virtual const inline void SetParameterFree() {
		fixed_ = false;
	}
	//! Set parameter free or fixed
	virtual const inline void FixParameter(const bool fixed) {
		fixed_ = fixed;
	}

	//! A public function returning a string naming its type
	/*!
	 * This function is used to get the type of the implementation of this
	 * general parameter interface.
	 * \sa operator<<, to_str(), make_str()
	 */
	virtual std::string TypeName() const {
		return "boolean";
	}

protected:
	//std::string out_; /*!< Output string to print information */
	bool hasError_; /*!< Is an error defined for this parameter? */
	bool usebounds_; /*!< Do you want to restrict your parameter? */
	bool fixed_; /*!< Do you want to keep parameter fixed? */
	int val_, err_; /*!< Containers of parameter information */

	//! A protected function which returns an output string for printing
	/*!
	 * This function uses all available information about the parameter
	 * to create a string which will be streamed via the stream operator <<.
	 * \sa operator<<, to_str(), type()
	 */
	virtual std::string make_str() {
		std::stringstream oss;
		oss << name_;
		oss << "\t Val = " << val_;
		if (hasError_)
			oss << "\t  Err = " << err_;
		oss << "\t Type = " << TypeName();
		return oss.str();
	}

	//! A protected function which returns an output string for printing
	/*!
	 * This function uses only the value information about the parameter
	 * to create a string which will be streamed via the stream operator <<.
	 * \sa operator<<, make_str()
	 */
	virtual std::string make_val_str() {
		std::stringstream ovs;
		ovs << val_;
		return ovs.str();
	}

};

} /* namespace ComPWA */

BOOST_SERIALIZATION_SHARED_PTR( ComPWA::DoubleParameter )
BOOST_CLASS_IMPLEMENTATION(
		ComPWA::DoubleParameter,
		boost::serialization::object_serializable
)
BOOST_CLASS_TRACKING(
		ComPWA::DoubleParameter,
		boost::serialization::track_never
)
BOOST_CLASS_IMPLEMENTATION(
		std::shared_ptr<ComPWA::DoubleParameter>,
		boost::serialization::object_serializable
)
BOOST_CLASS_TRACKING(
		std::shared_ptr<ComPWA::DoubleParameter>,
		boost::serialization::track_never
)
BOOST_CLASS_IMPLEMENTATION(
		std::vector<std::shared_ptr<ComPWA::DoubleParameter> >,
		boost::serialization::object_serializable
)
BOOST_CLASS_TRACKING(
		std::vector<std::shared_ptr<ComPWA::DoubleParameter> >,
		boost::serialization::track_never
)
#endif
