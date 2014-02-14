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
#ifndef ROOTEFFICIENCY_HPP_
#define ROOTEFFICIENCY_HPP_

#include <vector>
#include <memory>

#include <memory>

#include "TEfficiency.h"
#include "TH2.h"
#include "Physics/DPKinematics/DalitzKinematics.hpp"
#include "Core/Efficiency.hpp"
/**
 *  \class DalitzHistEfficiency
 *  \brief Efficiency provided by a histogram
 */
class DalitzHistEfficiency : public Efficiency {
protected:
	std::shared_ptr<TEfficiency> effHist;
public:
	//! returns efficiency for current datapoint
	virtual double evaluate(std::vector<double> x);
	virtual double evaluate(dataPoint& point);
	//! Construct DalitzHistEfficiency from TEfficiency object
	DalitzHistEfficiency(TEfficiency* eff);
	//! Construct DalitzHistEfficiency from two TH2 objects for passed and total events
	DalitzHistEfficiency(TH2D* passed, TH2D* total);
	DalitzHistEfficiency(const DalitzHistEfficiency&);
	~DalitzHistEfficiency(){};

};
/**
 *  \class DalitzAngleHistEfficiency
 *  \brief Uses also TEfficiency object, but the variables are one invariant
 *  mass and the corresponding helicity angle. This avoids binning effects near the phsp boundaries.
 *  ATTENTION: We assume that the invariant mass of particle 2 and 3 and the helicity angle
 *  between 1 and 2 are used!
 */
class DalitzAngleHistEfficiency : public DalitzHistEfficiency {
public:
	DalitzAngleHistEfficiency(TEfficiency* eff) : DalitzHistEfficiency(eff) {};
	DalitzAngleHistEfficiency(TH2D* passed, TH2D* total) : DalitzHistEfficiency(passed, total) {};
	DalitzAngleHistEfficiency(const DalitzAngleHistEfficiency& p) : DalitzHistEfficiency(p) {};
	~DalitzAngleHistEfficiency(){};
	virtual double evaluate(dataPoint& point);

};
/**
 *  \class DalitzPolyEfficiency
 *  \brief Efficiency provided by a polynomial
 */
class DalitzPolyEfficiency : public Efficiency {
private:

public:
	double evaluate() { return 1; };
	DalitzPolyEfficiency(){};
	~DalitzPolyEfficiency(){};

};

#endif /* ROOTEFFICIENCY_HPP_ */
