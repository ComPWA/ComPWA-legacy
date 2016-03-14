//-------------------------------------------------------------------------------
// Copyright (c) 2014 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff -
//-------------------------------------------------------------------------------
/*! \class SimpleAngleEfficiency
 * @file SimpleAngleEfficiency.hpp
 * Implementation of Efficiency interface. Efficiency over the Dalitz plot is stored in a
 * SimpleEfficiency object.
 */

#ifndef SIMPLEANGLEEFFICIENCY_HPP_
#define SIMPLEANGLEEFFICIENCY_HPP_
#include "Core/Efficiency.hpp"
#include "Physics/DPKinematics/SimpleEfficiency.hpp"
#include "TH1.h"
#include <vector>
#include <memory>

namespace ComPWA {
namespace Physics {
namespace DPKinematics {

class SimpleAngleEfficiency : public Efficiency {

protected:
	std::shared_ptr<SimpleEfficiency> effHist;
public:
	//! returns efficiency for current datapoint
	virtual double evaluate(std::vector<double> x);
	virtual double evaluate(dataPoint& point);
	//! Construct SimpleAngleEfficiency from TEfficiency object
	SimpleAngleEfficiency(SimpleEfficiency* eff);
	//! Construct SimpleAngleEfficiency from two TH2 objects for passed and total events
	SimpleAngleEfficiency(TH1* passed, TH1* total);
	SimpleAngleEfficiency(const SimpleAngleEfficiency&);
	~SimpleAngleEfficiency(){};

};

} /* namespace DPKinematics */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* SIMPLEANGLEEFFICIENCY_HPP_ */
