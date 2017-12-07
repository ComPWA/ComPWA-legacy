#include "Tools/ParameterTools.hpp"
#include "Physics/ParticleList.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/HelicityFormalism/IncoherentIntensity.hpp"
#include "Tools/Generate.hpp"
#include "Tools/RootGenerator.hpp"
#include "Tools/RootPlot.hpp"
#include "Tools/FitFractions.hpp"
#include "Core/Generator.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Kinematics.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
//#include <boost/python.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/log/trivial.hpp>
#include <boost/progress.hpp>

#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
//using namespace ComPWA;
using ComPWA::Physics::HelicityFormalism::HelicityKinematics;
using ComPWA::Physics::HelicityFormalism::IncoherentIntensity;
using ComPWA::Optimizer::Minuit2::MinuitResult;

PYBIND11_MAKE_OPAQUE(ComPWA::PartList);

py::array_t<double> result_values(std::shared_ptr<ComPWA::FitResult> fitRes){
  ComPWA::ParameterList resPar = fitRes->GetFinalParameters();
  //std::vector<double> ret;
 // for(unsigned int i=0; i<resPar.GetNParameter(); i++){
  //  ret.push_back(resPar.GetDoubleParameterValue(i));
 // }

  std::size_t size = resPar.GetNParameter();
  double *foo = new double[size];
  for (std::size_t i = 0; i < size; i++) {
      foo[i] = resPar.GetDoubleParameterValue(i);
  }

  // Create a Python object that will free the allocated
  // memory when destroyed:
  py::capsule free_when_done(foo, [](void *f) {
      double *foo = reinterpret_cast<double *>(f);
      //std::cerr << "Element [0] = " << foo[0] << "\n";
      //std::cerr << "freeing memory @ " << f << "\n";
      delete[] foo;
  });

  return py::array_t<double>(
      {size}, // shape
      {8}, // C-style contiguous strides for double
      foo, // the data pointer
      free_when_done); // numpy array references this parent
}

std::shared_ptr<ComPWA::AmpIntensity> createIntens(
		std::string modelStr,
		std::shared_ptr<ComPWA::PartList> partL,
		std::shared_ptr<ComPWA::Kinematics> kin,
		std::shared_ptr<ComPWA::DataReader::Data> phspSample){
  std::stringstream modelStream;
  modelStream << modelStr;
  boost::property_tree::ptree modelTree;
  boost::property_tree::xml_parser::read_xml(modelStream, modelTree);
  auto intens = IncoherentIntensity::Factory(
    partL, kin, modelTree.get_child("IncoherentIntensity"));
  auto phspPoints =
      std::make_shared<std::vector<ComPWA::dataPoint>>(phspSample->GetDataPoints(kin));
  intens->SetPhspSample(phspPoints, phspPoints);
  return intens;
}

std::vector<std::pair<std::string, std::string>> fitComponents(){
  std::vector<std::pair<std::string, std::string>> components;
  components.push_back(
    std::pair<std::string, std::string>("myAmp", "jpsiGammaPiPi"));
  components.push_back(
    std::pair<std::string, std::string>("f2(1270)", "jpsiGammaPiPi"));
  return components;
}

ComPWA::ParameterList calculateFitFractions(std::shared_ptr<ComPWA::Kinematics> kin,
		std::shared_ptr<ComPWA::AmpIntensity> intens, std::shared_ptr<ComPWA::DataReader::Data> phspSample){
  auto phspPoints =
		  std::make_shared<std::vector<ComPWA::dataPoint>>(phspSample->GetDataPoints(kin));
  return ComPWA::Tools::CalculateFitFractions(kin, intens, phspPoints, fitComponents());
}

void calcFractionError(ComPWA::ParameterList& fitPar, std::shared_ptr<ComPWA::FitResult> result,
		ComPWA::ParameterList& fitFracs, std::shared_ptr<ComPWA::AmpIntensity> intens, std::shared_ptr<ComPWA::Kinematics> kin,
		std::shared_ptr<ComPWA::DataReader::Data> phspSample, int nSets){
  auto resultM = std::dynamic_pointer_cast<MinuitResult>(result);
  auto phspPoints =
		  std::make_shared<std::vector<ComPWA::dataPoint>>(phspSample->GetDataPoints(kin));
  ComPWA::Tools::CalcFractionError(fitPar, resultM->GetCovarianceMatrix(), fitFracs, kin,
	                       intens, phspPoints, nSets, fitComponents());
}

void saveResults(std::string file, std::shared_ptr<ComPWA::FitResult> result){
  std::ofstream ofs(file);
  boost::archive::xml_oarchive oa(ofs);
  std::shared_ptr<MinuitResult> resultM = std::dynamic_pointer_cast<MinuitResult>(result);
  oa << BOOST_SERIALIZATION_NVP(resultM);
}

void saveModel(std::string file, std::shared_ptr<ComPWA::PartList> partL, ComPWA::ParameterList& fitPar, std::shared_ptr<ComPWA::AmpIntensity> intens){
  ComPWA::UpdateParticleList(partL, fitPar);
  boost::property_tree::ptree ptout;
  ptout.add_child("ParticleList", ComPWA::SaveParticles(partL));
  std::shared_ptr<IncoherentIntensity> intensI = std::dynamic_pointer_cast<IncoherentIntensity>(intens);
  ptout.add_child("IncoherentIntensity", IncoherentIntensity::Save(intensI));
  boost::property_tree::xml_parser::write_xml(file, ptout, std::locale());
}

class DataPoints {
public:
	DataPoints(std::shared_ptr<ComPWA::DataReader::Data> data, std::shared_ptr<ComPWA::Kinematics> kin) : nEvents(data->GetNEvents()), nVars(0) {
		std::vector<ComPWA::dataPoint> dataVec =  data->GetDataPoints(kin);
		nVars = dataVec[0].Size();
		rawEvtData = new double[nEvents*(nVars+2)]; //vars + weight + efficiency
        for (unsigned int i=0; i<data->GetNEvents(); i++){
          for (unsigned int j=0; j<dataVec[i].Size(); j++){
        	rawEvtData[nVars*i+j] = dataVec[i].GetPoint()[j];
          }
          rawEvtData[nVars*i+dataVec[i].Size()] = dataVec[i].GetWeight();
          rawEvtData[nVars*i+dataVec[i].Size()+1] = dataVec[i].GetEfficiency();
        }
    }
    double *getRawEvtData() { return rawEvtData; }
    size_t getNEvents() const { return nEvents; }
    size_t getNVars() const { return nVars; }
private:
    size_t nEvents;
    size_t nVars;
    double *rawEvtData;
};
 
//BOOST_PYTHON_MODULE(Dalitz_ext)
PYBIND11_MODULE(PyComPWA, m)
{
//    using namespace boost::python;
	m.doc() = "pybind11 DalitzFit plugin"; // optional module docstring

	// Global Functions
	m.def("ReadParticles", (void (*) (std::shared_ptr<ComPWA::PartList>, std::string)) &ComPWA::ReadParticles);
	m.def("GetDefaultParticles", []() { return defaultParticleList; } );
	m.def("GetInitialState", [](int id) { std::vector<ComPWA::pid> initialState = {id}; return initialState; } );
	m.def("GetFinalState", [](int idA, int idB, int idC) { std::vector<ComPWA::pid> finalState = {idA, idB, idC}; return finalState; } );
	m.def("GetIncoherentIntensity",  (std::shared_ptr<ComPWA::AmpIntensity>  (*)
			(std::string, std::shared_ptr<ComPWA::PartList>, std::shared_ptr<ComPWA::Kinematics>, std::shared_ptr<ComPWA::DataReader::Data>))
			&createIntens);
	m.def("setErrorOnParameterList", (void (*) (ComPWA::ParameterList&, double, bool)) &setErrorOnParameterList);
	m.def("Generate", (bool (*) (int, std::shared_ptr<ComPWA::Kinematics>,
    		std::shared_ptr<ComPWA::Generator>, std::shared_ptr<ComPWA::AmpIntensity>, std::shared_ptr<ComPWA::DataReader::Data>,
			std::shared_ptr<ComPWA::DataReader::Data>, std::shared_ptr<ComPWA::DataReader::Data>)) &ComPWA::Tools::Generate);
	        //phspTrue = std::shared_ptr<ComPWA::DataReader::Data>()
	m.def("GeneratePhsp", (bool (*) (int, std::shared_ptr<ComPWA::Generator>,
            std::shared_ptr<ComPWA::DataReader::Data>)) &ComPWA::Tools::GeneratePhsp);
	m.def("fitComponents",  (std::vector< std::pair<std::string, std::string> > (*) ()) &fitComponents);
	m.def("saveResults",  (void (*) (std::string, std::shared_ptr<ComPWA::FitResult>)) &saveResults);
	m.def("saveModel",  (void (*) (std::string, std::shared_ptr<ComPWA::PartList>, ComPWA::ParameterList&, std::shared_ptr<ComPWA::AmpIntensity>)) &saveModel);
	m.def("CalculateFitFractions", (ComPWA::ParameterList (*)
	        (std::shared_ptr<ComPWA::Kinematics>, std::shared_ptr<ComPWA::AmpIntensity>,
	        		std::shared_ptr<ComPWA::DataReader::Data> ))
			&calculateFitFractions);
	m.def("CalcFractionError", (void (*)
	        (ComPWA::ParameterList&, std::shared_ptr<ComPWA::FitResult>, ComPWA::ParameterList&,
	        		std::shared_ptr<ComPWA::AmpIntensity>, std::shared_ptr<ComPWA::Kinematics>,
	        		std::shared_ptr<ComPWA::DataReader::Data>, int))
			&calcFractionError);
	m.def("result_values",  (py::array_t<double> (*) (std::shared_ptr<ComPWA::FitResult>)) &result_values);

	// ComPWA Interfaces
	py::class_<ComPWA::Kinematics, std::shared_ptr<ComPWA::Kinematics> >(m, "Kinematics");
	py::class_<ComPWA::Generator, std::shared_ptr<ComPWA::Generator> >(m, "Generator");
	py::class_<ComPWA::IEstimator, std::shared_ptr<ComPWA::IEstimator> >(m, "IEstimator");
	py::class_<ComPWA::Optimizer::Optimizer, std::shared_ptr<ComPWA::Optimizer::Optimizer> >(m, "Optimizer");
	py::class_<ComPWA::FitResult, std::shared_ptr<ComPWA::FitResult> >(m, "FitResult")
	  .def("GetFinalParameters", &ComPWA::FitResult::GetFinalParameters)
	;
	py::class_<ComPWA::AbsParameter, std::shared_ptr<ComPWA::AbsParameter> >(m, "AbsParameter");

	// ComPWA Classes
	py::class_<ComPWA::AmpIntensity, std::shared_ptr<ComPWA::AmpIntensity> >(m, "AmpIntensity")
	  .def("GetParameters", &ComPWA::AmpIntensity::GetParameters)
	;

	py::class_<ComPWA::DataReader::Data, std::shared_ptr<ComPWA::DataReader::Data> >(m, "Data")
      .def(py::init<>())
	;

	py::class_<ComPWA::dataPoint>(m, "dataPoint")
      .def(py::init<>())
	;

	py::class_<DataPoints>(m, "DataPoints", py::buffer_protocol())
	   .def_buffer([](DataPoints &dp) -> py::buffer_info {
	        return py::buffer_info(
	            dp.getRawEvtData(),                               /* Pointer to buffer */
	            sizeof(double),                          /* Size of one scalar */
	            py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
	            2,                                      /* Number of dimensions */
	            { dp.getNEvents(), size_t(dp.getNVars()+2) },                 /* Buffer dimensions */
	            { sizeof(double) * dp.getNEvents(),             /* Strides (in bytes) for each index */
	              sizeof(double) }
	        );
	    })
		.def(py::init<std::shared_ptr<ComPWA::DataReader::Data>, std::shared_ptr<ComPWA::Kinematics> >())
	;

	py::class_<ComPWA::PartList, std::shared_ptr<ComPWA::PartList> >(m, "PartList")
      .def(py::init<>())
	;

	py::class_<ComPWA::DoubleParameter, ComPWA::AbsParameter, std::shared_ptr<ComPWA::DoubleParameter> >(m, "DoubleParameter")
      .def(py::init<>())
      .def(py::init<std::string, const double, const double >())
	;

	py::class_<ComPWA::ParameterList>(m, "ParameterList")
      .def(py::init<>())
	  .def("GetNParameter", &ComPWA::ParameterList::GetNParameter)
	  .def("SetParameterValue", (void (ComPWA::ParameterList::*) (const unsigned int, const double)) &ComPWA::ParameterList::SetParameterValue)
	  .def("SetParameterValue", (void (ComPWA::ParameterList::*) (std::string, const double)) &ComPWA::ParameterList::SetParameterValue)
	  .def("AddParameter", (void (ComPWA::ParameterList::*) (std::shared_ptr<ComPWA::DoubleParameter>)) &ComPWA::ParameterList::AddParameter)
	  //.def("GetDoubleParameter", (std::shared_ptr<ComPWA::DoubleParameter> (ComPWA::ParameterList::*) (const unsigned int)) &ComPWA::ParameterList::GetDoubleParameter)
	  //.def("GetDoubleParameter", (std::shared_ptr<ComPWA::DoubleParameter> (ComPWA::ParameterList::*) (std::string))  &ComPWA::ParameterList::GetDoubleParameter)
	  //.def("GetDoubleParameter", py::overload_cast<const unsigned int>(&ComPWA::ParameterList::GetDoubleParameter, py::const_))
	  //.def("GetDoubleParameter", py::overload_cast<const std::string>(&ComPWA::ParameterList::GetDoubleParameter, py::const_))
	  .def("ToString", &ComPWA::ParameterList::to_str)
	;

	py::class_<HelicityKinematics, ComPWA::Kinematics, std::shared_ptr<HelicityKinematics> >(m, "HelicityKinematics")
      .def(py::init<std::shared_ptr<ComPWA::PartList>, std::vector<ComPWA::pid>, std::vector<ComPWA::pid> >())
	;

	py::class_<ComPWA::Tools::RootGenerator, ComPWA::Generator, std::shared_ptr<ComPWA::Tools::RootGenerator> >(m, "RootGenerator")
      .def(py::init<std::shared_ptr<ComPWA::PartList>, std::shared_ptr<ComPWA::Kinematics> >())
	;

	py::class_<ComPWA::Tools::RootPlot, std::shared_ptr<ComPWA::Tools::RootPlot> >(m, "RootPlot")
      .def(py::init<std::shared_ptr<ComPWA::Kinematics> >())
	  .def("SetData", (void (ComPWA::Tools::RootPlot::*) (std::shared_ptr<ComPWA::DataReader::Data>)) &ComPWA::Tools::RootPlot::SetData)
	  .def("SetPhspSample", (void (ComPWA::Tools::RootPlot::*) (std::shared_ptr<ComPWA::DataReader::Data>)) &ComPWA::Tools::RootPlot::SetPhspSample)
	  .def("SetFitAmp", &ComPWA::Tools::RootPlot::SetFitAmp)
	  .def("AddComponent", &ComPWA::Tools::RootPlot::AddComponent)
	  .def("Write", &ComPWA::Tools::RootPlot::Write)
	;

	py::class_<ComPWA::Estimator::MinLogLH, ComPWA::IEstimator, std::shared_ptr<ComPWA::Estimator::MinLogLH> >(m, "MinLogLH")
      .def(py::init<std::shared_ptr<ComPWA::Kinematics>, std::shared_ptr<ComPWA::AmpIntensity>,
    		  std::shared_ptr<ComPWA::DataReader::Data>, std::shared_ptr<ComPWA::DataReader::Data>,
			  std::shared_ptr<ComPWA::DataReader::Data>, unsigned int, unsigned int >())
	  .def("UseFunctionTree", &ComPWA::Estimator::MinLogLH::UseFunctionTree)
	;

	py::class_<ComPWA::Optimizer::Minuit2::MinuitIF, ComPWA::Optimizer::Optimizer, std::shared_ptr<ComPWA::Optimizer::Minuit2::MinuitIF> >(m, "MinuitIF")
      .def(py::init<std::shared_ptr<ComPWA::IEstimator>, ComPWA::ParameterList&>())
	  .def("SetHesse", &ComPWA::Optimizer::Minuit2::MinuitIF::SetHesse)
	  .def("exec", &ComPWA::Optimizer::Minuit2::MinuitIF::exec)
	;

	py::class_<ComPWA::Optimizer::Minuit2::MinuitResult, ComPWA::FitResult, std::shared_ptr<ComPWA::Optimizer::Minuit2::MinuitResult> >(m, "MinuitResult")
	  .def("SetFitFractions", &ComPWA::Optimizer::Minuit2::MinuitResult::SetFitFractions)
	  .def("Print", &ComPWA::Optimizer::Minuit2::MinuitResult::Print, py::arg("opt") = "")
	;

}


/*USAGE
 * export PYTHONPATH=$PYTHONPATH:YOUR_COMPWA_BUILD_DIR/lib
 * python3
 * >>> from Dalitz_ext import *
 */
