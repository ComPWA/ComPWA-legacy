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

using namespace Gem::Geneva;

GenevaIF::GenevaIF(std::shared_ptr<ControlParameter> theData, std::string inConfigFileDir)
  : _myData(theData),configFileDir(inConfigFileDir),parallelizationMode(GO2_DEF_DEFAULPARALLELIZATIONMODE),
    serMode(Gem::Common::SERIALIZATIONMODE_BINARY),clientMode(false),ip("localhost"),port(0){

}

GenevaIF::~GenevaIF(){

}

void GenevaIF::setServerMode(){
  parallelizationMode = Gem::Geneva::EXECMODE_BROKERAGE;
  serMode = Gem::Common::SERIALIZATIONMODE_BINARY;
  clientMode = false;
}

void GenevaIF::setClientMode(std::string serverip, unsigned int serverport){
  parallelizationMode = Gem::Geneva::EXECMODE_BROKERAGE;
  serMode = Gem::Common::SERIALIZATIONMODE_BINARY;
  clientMode = true;
  ip = serverip;
  port = serverport;
}

std::shared_ptr<FitResult> GenevaIF::exec(ParameterList& par) {
	std::shared_ptr<GenevaResult> result(new GenevaResult());
	//Go2::init();
	//Go2 go(argc, argv, configFile);
	//Go2 go( clientMode, serMode, ip, port,
	 //   (configFileDir+"Go2.json"), parallelizationMode, GO2_DEF_DEFAULTVERBOSE);
	Go2 go( (configFileDir+"Go2.json"));

	//---------------------------------------------------------------------
	// Initialize a client, if requested

    if(go.clientMode()) {
	  std::cout << "Geneva Client waiting for action!" << std::endl;
	  go.clientRun();
	  return result;
	}

	//---------------------------------------------------------------------
	// Add individuals and algorithms and perform the actual optimization cycle

	//Provide Parameter in Geneva-Style
	unsigned int NPar = par.GetNDouble(); //just doubles up to now, TODO
	double val[NPar], min[NPar], max[NPar], err[NPar];
	std::string names[NPar];
	for(unsigned int i=0; i<NPar; i++){
	  std::shared_ptr<DoubleParameter> dpar = par.GetDoubleParameter(i);
	  names[i] = dpar->GetName();
	  val[i] = dpar->GetValue();
	  if(dpar->HasBounds()){
	    min[i] = dpar->GetMinValue();
	    max[i] = dpar->GetMaxValue();
	  }else{
	    max[i] = 1.79768e+307;//std::numeric_limits<double>::max();
	    min[i] = -1.79768e+307;//-1*max[i];
	  }
	  if(dpar->HasError())
	    err[i] = dpar->GetError()->GetError();
	  else
	    err[i] = val[i];
	}

	// Make an individual known to the optimizer
	boost::shared_ptr<GStartIndividual> p(new GStartIndividual(_myData, NPar, names, val, min, max, err));
	go.push_back(p);

	// Add an evolutionary algorithm to the Go2 class.
        GEvolutionaryAlgorithmFactory ea((configFileDir+"GEvolutionaryAlgorithm.json"), parallelizationMode);
	go & ea();

	// Perform the actual optimization
	boost::shared_ptr<GStartIndividual>
		bestIndividual_ptr = go.optimize<GStartIndividual>();

	// Do something with the best result

	// Terminate
//	double result= bestIndividual_ptr->getBestKnownFitness();
	//double finalValue = bestIndividual_ptr->getBestKnownFitness();
	result->setResult(bestIndividual_ptr);
	//int whattodowiththisidontknow =  go.finalize(); //Go2::finalize();

        //write Parameters back
	ParameterList resultPar;
	bestIndividual_ptr->getPar(resultPar);
	for(unsigned int i=0; i<par.GetNDouble(); i++){  //TODO: better way, no cast or check type
	  par.GetDoubleParameter(i)->SetValue(resultPar.GetDoubleParameter(i)->GetValue());
	}

	return result;
}
