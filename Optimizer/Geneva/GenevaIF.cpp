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

// Geneva header files go here
#include "geneva/Go2.hpp"
#include <geneva/GEvolutionaryAlgorithmFactory.hpp>

// ComPWA header files go here
#include "Optimizer/Geneva/GenevaIF.hpp"

// The individual that should be optimized
#include "Optimizer/Geneva/GStartIndividual.hpp"

using namespace Gem::Geneva;

GenevaIF::GenevaIF(std::shared_ptr<ControlParameter> theData, std::string inConfigFileDir)
  : _myData(theData),configFileDir(inConfigFileDir),parallelizationMode(GO2_DEF_DEFAULPARALLELIZATIONMODE),
    serMode(GO2_DEF_SERIALIZATIONMODE),clientMode(GO2_DEF_CLIENTMODE),ip(GO2_DEF_IP),port(GO2_DEF_PORT){

}

GenevaIF::~GenevaIF(){

}

void GenevaIF::setServerMode(){
  parallelizationMode = Gem::Geneva::PARMODE_BROKERAGE;
  serMode = Gem::Common::SERIALIZATIONMODE_BINARY;
  clientMode = false;
}

void GenevaIF::setClientMode(std::string serverip, unsigned int serverport){
  parallelizationMode = Gem::Geneva::PARMODE_BROKERAGE;
  serMode = Gem::Common::SERIALIZATIONMODE_BINARY;
  clientMode = true;
  ip = serverip;
  port = serverport;
}

const double GenevaIF::exec(ParameterList& par) {
	Go2::init();
	//Go2 go(argc, argv, configFile);
	Go2 go( clientMode, serMode, ip, port,
	    (configFileDir+"Go2.json"), parallelizationMode, GO2_DEF_DEFAULTVERBOSE);

	//---------------------------------------------------------------------
	// Initialize a client, if requested

        if(go.clientMode()) {
	  std::cout << "Geneva Client waiting for action!" << std::endl;
	  return go.clientRun();
	}

	//---------------------------------------------------------------------
	// Add individuals and algorithms and perform the actual optimization cycle

	//Provide Parameter in Geneva-Style
	unsigned int NPar = par.GetNDouble(); //just doubles up to now, TODO
	double val[NPar], min[NPar], max[NPar], err[NPar];
	std::string names[NPar];
	for(unsigned int i=0; i<NPar; i++){
	  names[i] = par.GetDoubleParameter(i).GetName();
	  val[i] = par.GetDoubleParameter(i).GetValue();
	  min[i] = par.GetDoubleParameter(i).GetMinValue();
	  max[i] = par.GetDoubleParameter(i).GetMaxValue();
	  err[i] = par.GetDoubleParameter(i).GetError();
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
	double result = bestIndividual_ptr->getBestKnownFitness();
	int whattodowiththisidontknow =  go.finalize(); //Go2::finalize();

        //write Parameters back
	ParameterList resultPar;
	bestIndividual_ptr->getPar(resultPar);
	for(unsigned int i=0; i<par.GetNDouble(); i++){  //TODO: better way, no cast or check type
	  par.GetDoubleParameter(i).SetValue(resultPar.GetDoubleParameter(i).GetValue());
	}

	return result;
}
