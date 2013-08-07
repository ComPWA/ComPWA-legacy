/*Author: Peter Weidenkaff
 * Date: 2013-08-07
 */
#ifndef PLOTDATA_HPP_
#define PLOTDATA_HPP_

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TObject.h>
#include <TCanvas.h>

//Core header files go here
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"

// ComPWA header files go here
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"

class plotData
{
public:
	//	plotData();
	plotData(std::string,std::string,std::shared_ptr<Data>,std::shared_ptr<Data>,std::shared_ptr<Amplitude>);
	plotData(std::string,std::string,std::shared_ptr<Data>,std::shared_ptr<Data>);
	plotData(std::string,std::string,std::shared_ptr<Data>);
	void setStyle(unsigned int);
	void plot();
	void setPar(ParameterList par){_par=par;};
	//  virtual ~plotData();
protected:
	std::shared_ptr<Data> _set1;
	std::shared_ptr<Data> _set2;
	std::shared_ptr<Amplitude> _amp;
	std::string _outFile;
	ParameterList _par;
	std::string _name;
	unsigned int _style;

};

#endif
