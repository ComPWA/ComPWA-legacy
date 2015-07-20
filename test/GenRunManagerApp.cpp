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

#include "../Physics/DPKinematics/RootGenerator.hpp"
// Physics Interface header files go here
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/RunManager.hpp"
#include "Core/Efficiency.hpp"
#include "DataReader/RootReader/RootReader.hpp"

#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"

//#include "PWA/PlotData.hpp"

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(DalitzKinematics::createInstance("J/psi","gamma","pi0","pi0"));
	std::string outFile="3Part-4vecs.root";
	unsigned int dataSize = 100000;
	//load resonances
	AmplitudeSetup ini("test/JPSI_ypipi.xml");
	cout << "loaded file " << ini.getFileName() << " with " << ini.getBreitWigner().size() << " resonances!" << std::endl;
	std::shared_ptr<Data> data(new RootReader(outFile, "data",true));
	std::shared_ptr<Data> phsp(new RootReader(outFile, "mc",true));
	std::shared_ptr<Generator> gen(new RootGenerator());
	std::shared_ptr<Amplitude> amp(new AmpSumIntensity(
			ini,std::shared_ptr<Efficiency>(new UnitEfficiency()),normStyle::none));

	RunManager run(dataSize,amp,gen);
	run.setGenerator(gen);
	run.setData(data);
	run.setPhspSample(phsp);
	run.generate(dataSize);
	run.generatePhsp(dataSize);
	std::cout<<"Data size: "<<data->getNEvents()<<std::endl;
	data->writeData();
	phsp->writeData();
//	TFile output(outFile.c_str(),"update");
//	output.SetCompressionLevel(1); //try level 2 also
//	output.Close(); //clean output file
//	plotData plot(std::string("dalitz"),outFile, data);  plot.plot(300);
//	plotData plotPhsp(std::string("dalitz_phsp"),outFile, data);  plotPhsp.plot(300);

	return 0;
}
