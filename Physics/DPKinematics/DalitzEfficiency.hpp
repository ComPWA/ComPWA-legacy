//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//		Peter Weidenkaff -
//-------------------------------------------------------------------------------
#ifndef DALITZEFFICIENCY_HPP_
#define DALITZEFFICIENCY_HPP_

#include "TEfficiency.h"
#include "TH2.h"
/**
 *  \class DalitzEfficiency
 *  \brief Virtual efficiency class
 */
class DalitzEfficiency{
private:

public:
  DalitzEfficiency();

  DalitzEfficiency(const DalitzEfficiency&, const char*);

  virtual ~DalitzEfficiency();

	virtual double evaluate() = 0;
};

/**
 *  \class DalitzHistEfficiency
 *  \brief Efficiency provided by a histogram
 */
class DalitzHistEfficiency : public DalitzEfficiency {
private:
	std::shared_ptr<TEfficiency> effHist;
public:
	//! returns efficiency for current datapoint
	double evaluate();
	//! Construct DalitzHistEfficiency from TEfficiency object
	DalitzHistEfficiency(TEfficiency* eff);
	//! Construct DalitzHistEfficiency from two TH2 objects for passed and total events
	DalitzHistEfficiency(TH2D* passed, TH2D* total);
	DalitzHistEfficiency(const DalitzEfficiency&);
	~DalitzHistEfficiency(){};

};
/**
 *  \class DalitzPolyEfficiency
 *  \brief Efficiency provided by a polynomial
 */
class DalitzPolyEfficiency : public DalitzEfficiency {
private:

public:
	double evaluate() { return 1; };
	DalitzPolyEfficiency(){};
	~DalitzPolyEfficiency(){};

};

#endif /* DALITZEFFICIENCY_HPP_ */
