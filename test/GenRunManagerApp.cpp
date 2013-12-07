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

// Physics Interface header files go here
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "DataReader/RootReader/RootReader.hpp"

#include "Physics/DPKinematics/DalitzKinematics.hpp"
#include "Core/RunManager.hpp"
#include "Physics/DPKinematics/DataPoint.hpp"
#include "Physics/DPKinematics/RootGenerator.cpp"
#include "Core/Efficiency.hpp"

//#include "PWA/PlotData.hpp"

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
//	static dataPoint* point = dataPoint::instance();
	DalitzKinematics* kin = DalitzKinematics::createInstance("J/psi","gamma","pi0","pi0");
	std::string outFile="gen.root";
	unsigned int dataSize = 10000;
	//load resonances
	AmplitudeSetup ini("test/JPSI_ypipi.xml");
	cout << "loaded file " << ini.getFileName() << " with " << ini.getResonances().size() << " resonances!" << std::endl;
	std::shared_ptr<Data> data(new RootReader(outFile, true,"lalelu",false));
	std::shared_ptr<Generator> gen(new RootGenerator());
	std::shared_ptr<Amplitude> amp(new AmpSumIntensity(ini,dataSize,AmpSumIntensity::normalizationStyle::one));
	std::shared_ptr<Efficiency> eff(new UnitEfficiency());


	RunManager run(dataSize,amp,eff,gen);
	run.setData(data);
	run.generate(dataSize);
	std::cout<<"Data size: "<<data->getNEvents()<<std::endl;
//	TFile output(outFile.c_str(),"update");
//	output.SetCompressionLevel(1); //try level 2 also
//	output.Close(); //clean output file
//	plotData plot(std::string("dalitz"),outFile, data);  plot.plot(300);
//	plotData plotPhsp(std::string("dalitz_phsp"),outFile, data);  plotPhsp.plot(300);

	return 0;
}
