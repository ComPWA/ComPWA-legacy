//! Application to generate J/Psi -> g pi pi.
/*!
 * @file GenDalitzApp.cpp
 * This application uses the Breit-Wigner-Amplitude-Sum module and a
 * phase-space generator to generate a file with J/Psi -> g pi pi events.
 */

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

<<<<<<< HEAD

#include "../Physics/DPKinematics/RootGenerator.hpp"
=======
>>>>>>> dd7eb340b0f73d2b07005f687e586aae14fbc9fa
// Physics Interface header files go here
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/RunManager.hpp"
#include "Core/Efficiency.hpp"
#include "DataReader/RootReader/RootReader.hpp"
#include "DataReader/RootGenerator/RootGenerator.hpp"

#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"

//#include "PWA/PlotData.hpp"

using namespace ComPWA;
using Physics::DPKinematics::DalitzKinematics;
using DataReader::RootReader::RootReader;
using DataReader::RootGenerator::RootGenerator;
using DataReader::Data;
using Physics::AmplitudeSum::AmplitudeSetup;
using Physics::AmplitudeSum::AmpSumIntensity;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(DalitzKinematics::createInstance("J/psi","gamma","pi0","pi0"));

	std::string outFile="3Part-4vecs.root";
	unsigned int dataSize = 100000;

	std::shared_ptr<Data> data(new RootReader(outFile, "data",true));
	std::shared_ptr<Data> phsp(new RootReader(outFile, "mc",true));
	std::shared_ptr<Generator> gen(new RootGenerator());

	std::string resoFile = "/test/JPSI_ypipi.xml";
	boost::property_tree::ptree pt;
	read_xml(resoFile , pt, boost::property_tree::xml_parser::trim_whitespace);
	auto a = new AmpSumIntensity("amp",normStyle::none,
			std::shared_ptr<Efficiency>(new UnitEfficiency()), dataSize);
	a->Configure(pt);
	a->to_str();
	std::shared_ptr<Amplitude> amp(a);

	RunManager run(dataSize,amp,gen);
	run.setGenerator(gen);
	run.setData(data);
	run.setPhspSample(phsp);
	run.generate(dataSize);
	run.generatePhsp(dataSize);
	std::cout<<"Data size: "<<data->getNEvents()<<std::endl;
	data->writeData();
	phsp->writeData();

/* ****Pull generation by Mathias, TODO: better handling
	for(unsigned int i=0; i<1; i++){
	  std::string outFile="test/3Part-4vecs_1M_PDG_PHASE_";
	  outFile+=std::to_string(i);
	  outFile+=".root";
	  unsigned int dataSize = 1000000;
	  //load resonances
	  AmplitudeSetup ini("test/JPSI_ypipi_phase.xml");
	  cout << "loaded file " << ini.getFileName() << " with " << ini.getResonances().size() << " resonances!" << std::endl;
	  std::shared_ptr<Data> data(new RootReader(outFile, true,"data",false));
	  std::shared_ptr<Data> phsp(new RootReader(outFile, true,"mc",false));
	  std::shared_ptr<Generator> gen(new RootGenerator(i*7+13));
	  std::shared_ptr<Amplitude> amp(new AmpSumIntensity(ini,std::shared_ptr<Efficiency>(new UnitEfficiency()),dataSize,AmpSumIntensity::normStyle::none));

	  RunManager run(dataSize,amp,gen);
	  run.setGenerator(gen);
	  run.setData(data);
	  run.setPhspSample(phsp);
	  run.generate(dataSize);
	  run.generatePhsp(dataSize);
	  std::cout<<"Data size: "<<data->getNEvents()<<std::endl;
	  data->writeData();
	  phsp->writeData();
	}*/

//	TFile output(outFile.c_str(),"update");
//	output.SetCompressionLevel(1); //try level 2 also
//	output.Close(); //clean output file
//	plotData plot(std::string("dalitz"),outFile, data);  plot.plot(300);
//	plotData plotPhsp(std::string("dalitz_phsp"),outFile, data);  plotPhsp.plot(300);

	return 0;
}
