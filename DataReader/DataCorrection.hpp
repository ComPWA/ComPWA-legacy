/*
 * DataCorrection.hpp
 *
 *  Created on: Aug 27, 2015
 *      Author: weidenka
 */

#ifndef DATAREADER_DATACORRECTION_HPP_
#define DATAREADER_DATACORRECTION_HPP_

#include "Core/Event.hpp"
#include "Core/Logging.hpp"
#include "DataReader/CorrectionTable.hpp"
#include <stdexcept>
#include <cfloat>

namespace ComPWA {


class DataCorrection
{
public:
	virtual ~DataCorrection() { }
	virtual double getCorrection(Event& ev) = 0;
	virtual void Print() const = 0;
};

class UnitCorrection : public DataCorrection
{
	virtual ~UnitCorrection() { }
	virtual double getCorrection(Event& ev) { return 1; }
	virtual void Print() const {
		BOOST_LOG_TRIVIAL(info) << "UnitCorrection::Print() | correction factors are set to one";
	}
};

class MomentumCorrection : public DataCorrection
{
public:
	MomentumCorrection(std::vector<CorrectionTable> inCorr, std::string t="");
	~MomentumCorrection() { };
	double getCorrection(Event& ev);
	virtual void Print() const;
	//! Get title
	std::string GetTitle() { return title; }
	//! Set title
	void SetTitle( std::string t ) { title=t; }
protected:
	std::vector<CorrectionTable> corrections;
	std::string title;
};

}
#endif /* DATAREADER_DATACORRECTION_HPP_ */
