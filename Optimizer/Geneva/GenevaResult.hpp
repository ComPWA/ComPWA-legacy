//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//			Peter Weidenkaff
//          Mathias Michel
//-------------------------------------------------------------------------------
//! Optimizer Interface Base-Class.
/*! \class GenevaResult
 * @file GenevaResult.hpp
 * This class contains FitResults from Geneva Optimizations
 */

#ifndef _GENEVARESULT_HPP_
#define _GENEVARESULT_HPP_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/shared_ptr.hpp>

#include "Core/ParameterList.hpp"
#include "Core/TableFormater.hpp"
#include "Core/PhysConst.hpp"
#include "Core/FitResult.hpp"

#include "Optimizer/Geneva/GStartIndividual.hpp"

class GenevaResult : public FitResult
{
public:
	GenevaResult() : finalLH(0) {};
	void setResult(boost::shared_ptr<Gem::Geneva::GStartIndividual> result){ init(result); }
	operator double() const { return finalLH; };
	double getResult(){ return finalLH; }

protected:
	double finalLH;

	void genOutput(std::ostream& out, std::string opt="");
	void init(boost::shared_ptr<Gem::Geneva::GStartIndividual> min);
};

#endif
