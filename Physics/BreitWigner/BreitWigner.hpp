//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Physics Module with simple 1D Breit-Wigner.
/*! \class BreitWigner
 * @file BreitWigner.hpp
 * This class provides a simple Breit-Wigner calculation with given parameters.
 * It fulfills the Amplitude interface to be compatible with other ComPWA modules.
 */

#ifndef _PIFBW_HPP
#define _PIFBW_HPP

#include <vector>
#include <memory>
#include "Core/Amplitude.hpp"
#include "Core/ParameterList.hpp"

class BreitWigner : public Amplitude {

public:
	/// Default Constructor (0x0)
	BreitWigner(const double min, const double max);

	//! normalization integral for parameters \par
	virtual const double normalization(ParameterList& par) { return integral(par); }
	//! normalization integral
	virtual const double normalization() { return integral(); };
	//! integral for parameters \par excluding efficiency
	virtual const double integral() { };
	//! integral excluding efficiency
	virtual const double integral(ParameterList& par);

	virtual double getMaxVal(ParameterList& par, std::shared_ptr<Generator> gen ){ return 1; };
	virtual double getMaxVal(std::shared_ptr<Generator> gen) { return 1; };

	virtual const double volume();
	virtual const double drawInt(double *x, double *p); //For easy usage in a root TF1
	virtual const ParameterList& intensity(double x, double M, double T);
	virtual const ParameterList& intensity(std::vector<double> x, ParameterList& par);
	virtual const ParameterList& intensity(dataPoint& point, ParameterList& par);
	virtual const ParameterList& intensity(dataPoint& point);
	virtual const ParameterList& intensityNoEff(dataPoint& point){ return intensity(point); };
	virtual bool copyParameterList(ParameterList& outPar);

	virtual void setParameterList(ParameterList& par);
	virtual void setNevents(unsigned int n) { _entries=n; };
	virtual unsigned int getNevents() { return _entries; };
	/** Destructor */
	virtual ~BreitWigner();
	virtual BreitWigner* Clone() {
		return new BreitWigner(*this);
	}

	/* not implemented yet */
	virtual void printAmps() {};
	virtual double getFraction(unsigned int id) { return 1.;};
	virtual double getFraction(std::string name) { return 1.;};
	virtual void printFractions() {};
	virtual double getIntValue(std::string var1, double min1, double max1, std::string var2, double min2, double max2) {};

protected:
	double min_;
	double max_;
	unsigned int _entries;
	ParameterList params;

private:
	const double BreitWignerValue(double x, double M, double T);

};

#endif /* _PIFBW_HPP */
