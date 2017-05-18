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
/*
 * Copyright (C) Gemfony scientific UG (haftungsbeschraenkt)
 *
 * See the AUTHORS file in the top-level directory for a list of authors.
 *
 * Contact: contact [at] gemfony (dot) com
 *
 * This file is part of the Geneva library collection.
 *
 * Geneva was developed with kind support from Karlsruhe Institute of
 * Technology (KIT) and Steinbuch Centre for Computing (SCC). Further
 * information about KIT and SCC can be found at http://www.kit.edu/english
 * and http://scc.kit.edu .
 *
 * Geneva is free software: you can redistribute and/or modify it under
 * the terms of version 3 of the GNU Affero General Public License
 * as published by the Free Software Foundation.
 *
 * Geneva is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with the Geneva library. If not, see <http://www.gnu.org/licenses/>.
 *
 * For further information on Gemfony scientific and Geneva, visit
 * http://www.gemfony.com .
 */

// Standard header files go here
#include <iostream>
#include <limits>
#include <memory>

// Geneva header files go here
#include "geneva/Go2.hpp"
#include <geneva/GEvolutionaryAlgorithmFactory.hpp>

// ComPWA header files go here
#include "Optimizer/Geneva/GenevaIF.hpp"
#include "Optimizer/Geneva/GenevaResult.hpp"

// The individual that should be optimized
#include "Optimizer/Geneva/GStartIndividual.hpp"

namespace ComPWA {
namespace Optimizer {
namespace Geneva {

using namespace Gem::Geneva;

GenevaIF::GenevaIF(std::shared_ptr<ControlParameter> theData, std::string inConfigFileDir)
  : _myData(theData),configFileDir(inConfigFileDir),parallelizationMode(execMode::EXECMODE_SERIAL),
    serMode(Gem::Common::serializationMode::SERIALIZATIONMODE_BINARY),clientMode(false),ip("localhost"),port(0){
	BOOST_LOG_TRIVIAL(info) << "GenevaIF::GenevaIF() | "
			"Starting Geneva interface: config dir="<<configFileDir;
}

GenevaIF::~GenevaIF(){

}

void GenevaIF::setServerMode(){
  parallelizationMode = execMode::EXECMODE_BROKERAGE;
  serMode = Gem::Common::serializationMode::SERIALIZATIONMODE_BINARY;
  clientMode = false;
}

void GenevaIF::setClientMode(std::string serverip, unsigned int serverport){
  parallelizationMode = execMode::EXECMODE_BROKERAGE;
  serMode = Gem::Common::serializationMode::SERIALIZATIONMODE_BINARY;
  clientMode = true;
  ip = serverip;
  port = serverport;
}

std::shared_ptr<FitResult> GenevaIF::exec(ParameterList& par) {
	std::shared_ptr<GenevaResult> result(new GenevaResult());
	ParameterList initialParList(par);

	Go2 go( (configFileDir+"Go2.json"));

	// Initialize a client, if requested
    if(go.clientMode()) {
	  std::cout << "Geneva Client waiting for action!" << std::endl;
	  go.clientRun();
	  return result;
	}

	std::shared_ptr<GStartIndividual> p( new GStartIndividual(_myData, par) );
	go.push_back(p);

	// Add an evolutionary algorithm to the Go2 class.
	GEvolutionaryAlgorithmFactory ea((configFileDir+"GEvolutionaryAlgorithm.json"), parallelizationMode);
	go & ea();

	// Perform the actual optimization
	std::shared_ptr<GStartIndividual>
		bestIndividual_ptr = go.optimize<GStartIndividual>();

	// Terminate
	bestIndividual_ptr->getPar(par);
	result->setResult(bestIndividual_ptr);
	//result->SetAmplitude(_myData->getAmplitudes().at(0));
	//result->setInitialParameters(initialParList);
	//result->setFinalParameters(par);
	//int whattodowiththisidontknow =  go.finalize(); //Go2::finalize();

        //write Parameters back
	ParameterList resultPar;
	bestIndividual_ptr->getPar(resultPar);
	for(unsigned int i=0; i<par.GetNDouble(); i++){  //TODO: better way, no cast or check type
	  if(!par.GetDoubleParameter(i)->IsFixed())
	    par.GetDoubleParameter(i)->SetValue(resultPar.GetDoubleParameter(i)->GetValue());
	}

	return result;
}

} /* namespace Geneva */
} /* namespace Optimizer */
} /* namespace ComPWA */
