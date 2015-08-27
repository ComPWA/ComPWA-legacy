/*
 * DataCorrection.hpp
 *
 *  Created on: Aug 27, 2015
 *      Author: weidenka
 */

#ifndef DATAREADER_DATACORRECTION_HPP_
#define DATAREADER_DATACORRECTION_HPP_

#include "Core/Event.hpp"
#include "DataReader/CorrectionTable.hpp"
#include <stdexcept>
#include <cfloat>

class DataCorrection
{
public:
	virtual ~DataCorrection() { }
	virtual double getCorrection(Event& ev) = 0;
};

class UnitCorrection : public DataCorrection
{
	virtual ~UnitCorrection() { }
	virtual double getCorrection(Event& ev) { return 1; }
};

class MomentumCorrection : public DataCorrection
{
public:
	MomentumCorrection(std::vector<CorrectionTable> inCorr);
	~MomentumCorrection() { };
	double getCorrection(Event& ev);
protected:
	std::vector<CorrectionTable> corrections;
};



#endif /* DATAREADER_DATACORRECTION_HPP_ */
