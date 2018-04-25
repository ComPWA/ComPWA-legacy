#include "Tools/ParameterTools.hpp"
#include "Physics/ParticleList.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IncoherentIntensity.hpp"
#include "Tools/RootGenerator.hpp"
#include "Tools/RootPlot.hpp"
#include "Tools/FitFractions.hpp"
#include "Tools/Generate.hpp"
#include "Core/Generator.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Kinematics.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

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
// using namespace ComPWA;
using ComPWA::Physics::HelicityFormalism::HelicityKinematics;
using ComPWA::Physics::IncoherentIntensity;
using ComPWA::Optimizer::Minuit2::MinuitResult;

PYBIND11_MAKE_OPAQUE(ComPWA::PartList);
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

/// Helper function to create a Incoherent intensity from a model file.
std::shared_ptr<ComPWA::AmpIntensity>
createIntens(std::string modelStr, std::shared_ptr<ComPWA::PartList> partL,
             std::shared_ptr<ComPWA::Kinematics> kin,
             std::shared_ptr<ComPWA::DataReader::Data> phspSample,
             std::shared_ptr<ComPWA::DataReader::Data> truePhspSample) {
  std::stringstream modelStream;
  modelStream << modelStr;
  boost::property_tree::ptree modelTree;
  boost::property_tree::xml_parser::read_xml(modelStream, modelTree);
  auto intens = std::make_shared<IncoherentIntensity>(
      partL, kin, modelTree.get_child("Intensity"));

  // Setting phsp samples. The true sample does not include detector efficiency
  // and is needed to calculate the normalization of components.
  if (phspSample == truePhspSample) {
    auto phspPoints = std::make_shared<std::vector<ComPWA::DataPoint>>(
        phspSample->dataPoints(kin));
    intens->setPhspSample(phspPoints, phspPoints);
  } else {
    auto truePhspPoints = std::make_shared<std::vector<ComPWA::DataPoint>>(
        truePhspSample->dataPoints(kin));
    auto phspPoints = std::make_shared<std::vector<ComPWA::DataPoint>>(
        phspSample->dataPoints(kin));
    intens->setPhspSample(phspPoints, truePhspPoints);
  }
  return intens;
}

class DataPoints {
public:
  DataPoints(std::shared_ptr<ComPWA::DataReader::Data> data,
             std::shared_ptr<ComPWA::Kinematics> kin)
      : nEvents(data->numEvents()), nVars(0) {
    std::vector<ComPWA::DataPoint> dataVec = data->dataPoints(kin);
    nVars = dataVec[0].size();
    rawEvtData = new double[nEvents * (nVars + 2)]; // vars + weight +
                                                    // efficiency
    for (unsigned int i = 0; i < data->numEvents(); i++) {
      for (unsigned int j = 0; j < dataVec[i].size(); j++) {
        rawEvtData[nVars * i + j] = dataVec[i].values()[j];
      }
      rawEvtData[nVars * i + dataVec[i].size()] = dataVec[i].weight();
      rawEvtData[nVars * i + dataVec[i].size() + 1] = dataVec[i].efficiency();
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

///
/// Python interface for ComPWA
///
/// USAGE:
/// export PYTHONPATH=$PYTHONPATH:YOUR_COMPWA_BUILD_DIR/lib
/// python3
/// >>> from PyComPWA import *
///
PYBIND11_MODULE(pycompwa, m) {
  m.doc() = "ComPWA python interface"; // optional module docstring

  // ComPWA interface classes
  py::class_<ComPWA::Kinematics, std::shared_ptr<ComPWA::Kinematics>>(
      m, "Kinematics");

  py::class_<ComPWA::Generator, std::shared_ptr<ComPWA::Generator>>(
      m, "Generator");

  py::class_<ComPWA::IEstimator, std::shared_ptr<ComPWA::IEstimator>>(
      m, "Estimator");

  py::class_<ComPWA::Optimizer::Optimizer,
             std::shared_ptr<ComPWA::Optimizer::Optimizer>>(m, "Optimizer");

  py::class_<ComPWA::FitResult, std::shared_ptr<ComPWA::FitResult>>(m,
                                                                    "FitResult")
      .def("final_parameters", &ComPWA::FitResult::finalParameters);

  py::class_<ComPWA::Parameter, std::shared_ptr<ComPWA::Parameter>>(
      m, "Parameter");

  py::class_<ComPWA::DataReader::Data,
             std::shared_ptr<ComPWA::DataReader::Data>>(m, "Data")
      .def(py::init<>())
      .def("size", &ComPWA::DataReader::Data::numEvents)
      .def("num_events", &ComPWA::DataReader::Data::numEvents);

  py::class_<ComPWA::DataPoint>(m, "DataPoint").def(py::init<>());

  // Global Functions
  m.def("log", [](std::string msg) { LOG(info) << msg; });

  m.def("read_particles",
        (void (*)(std::shared_ptr<ComPWA::PartList>, std::string)) &
            ComPWA::ReadParticles);

  m.def("default_particles", []() { return defaultParticleList; });

  m.def("incoherent_intensity",
        (std::shared_ptr<ComPWA::AmpIntensity>(*)(
            std::string, std::shared_ptr<ComPWA::PartList>,
            std::shared_ptr<ComPWA::Kinematics>,
            std::shared_ptr<ComPWA::DataReader::Data>,
            std::shared_ptr<ComPWA::DataReader::Data>)) &
            createIntens);

  m.def("set_parameter_error",
        (void (*)(ComPWA::ParameterList &, double, bool)) &
            setErrorOnParameterList);

  m.def("generate", (bool (*)(int, std::shared_ptr<ComPWA::Kinematics>,
                              std::shared_ptr<ComPWA::Generator>,
                              std::shared_ptr<ComPWA::AmpIntensity>,
                              std::shared_ptr<ComPWA::DataReader::Data>,
                              std::shared_ptr<ComPWA::DataReader::Data>,
                              std::shared_ptr<ComPWA::DataReader::Data>)) &
                        ComPWA::Tools::generate,
        "Generate sample from AmpIntensity", py::arg("size"), py::arg("kin"),
        py::arg("gen"), py::arg("intens"), py::arg("sample"),
        py::arg("phspSample") = std::shared_ptr<ComPWA::DataReader::Data>(),
        py::arg("toyPhspSample") = std::shared_ptr<ComPWA::DataReader::Data>());

  m.def("generate_phsp", (bool (*)(int, std::shared_ptr<ComPWA::Generator>,
                                   std::shared_ptr<ComPWA::DataReader::Data>)) &
                             ComPWA::Tools::generatePhsp);

  m.def("fit_fractions",
        [](std::shared_ptr<ComPWA::Kinematics> kin,
           std::shared_ptr<ComPWA::AmpIntensity> intens,
           std::shared_ptr<ComPWA::DataReader::Data> toyPhspSample,
           std::vector<std::pair<std::string, std::string>> components) {
          auto toyPhspPoints = std::make_shared<std::vector<ComPWA::DataPoint>>(
              toyPhspSample->dataPoints(kin));
          return ComPWA::Tools::CalculateFitFractions(
              kin, intens, toyPhspPoints, components);
        },
        "Calculate fit fractions for a list of components given an "
        "intensity.",
        py::arg("kin"), py::arg("intensity"), py::arg("sample"),
        py::arg("components"));
  
  m.def(
      "fit_fractions_error",
      [](ComPWA::ParameterList &fitParameters,
         std::shared_ptr<ComPWA::FitResult> fitResult,
         ComPWA::ParameterList &fitFractions,
         std::shared_ptr<ComPWA::AmpIntensity> intens,
         std::vector<std::pair<std::string, std::string>> components,
         std::shared_ptr<ComPWA::Kinematics> kin,
         std::shared_ptr<ComPWA::DataReader::Data> phspSample, int nSets) {
        auto resultM = std::dynamic_pointer_cast<MinuitResult>(fitResult);
        auto phspPoints = std::make_shared<std::vector<ComPWA::DataPoint>>(
            phspSample->dataPoints(kin));
        ComPWA::Tools::CalcFractionError(
            fitParameters, resultM->covarianceMatrix(), fitFractions, kin,
            intens, phspPoints, nSets, components);
      },
      "Calculate uncertainties fot a list of fit fractions given a fit result "
      "and its corresponding intensity.",
      py::arg("fit_params"), py::arg("fit_result"), py::arg("fit_fractions"),
      py::arg("intensity"), py::arg("components"), py::arg("kin"),
      py::arg("phspSample"), py::arg("nSets"));

  // ComPWA Classes
  py::class_<ComPWA::Logging, std::shared_ptr<ComPWA::Logging>>(m, "Logging")
      .def(py::init<std::string, std::string>());

  py::class_<ComPWA::AmpIntensity, std::shared_ptr<ComPWA::AmpIntensity>>(
      m, "AmpIntensity")
      .def("parameters", &ComPWA::AmpIntensity::parameters)
      .def("write",
           [](const ComPWA::AmpIntensity &intens, std::string file) {
             boost::property_tree::ptree ptout;
             ptout.add_child("IncoherentIntensity", intens.save());
             boost::property_tree::xml_parser::write_xml(file, ptout,
                                                         std::locale());
           },
           py::arg("file"));

  py::class_<DataPoints>(m, "DataPoints", py::buffer_protocol())
      .def_buffer([](DataPoints &dp) -> py::buffer_info {
        return py::buffer_info(
            dp.getRawEvtData(),                      // Pointer to buffer
            sizeof(double),                          // Size of one scalar
            py::format_descriptor<double>::format(), // Python struct-style
                                                     // format descriptor
            2,                                       // Number of dimensions
            {dp.getNEvents(), size_t(dp.getNVars() + 2)}, // Buffer dimensions
            {sizeof(double) *
                 dp.getNEvents(), // Strides (in bytes) for each index
             sizeof(double)});
      })
      .def(py::init<std::shared_ptr<ComPWA::DataReader::Data>,
                    std::shared_ptr<ComPWA::Kinematics>>());

  py::class_<ComPWA::PartList, std::shared_ptr<ComPWA::PartList>>(m, "PartList")
      .def(py::init<>())
      .def("write",
           [](const ComPWA::PartList &list, std::string file) {
             boost::property_tree::ptree ptout;
             ptout.add_child("ParticleList", SaveParticles(list));
             boost::property_tree::xml_parser::write_xml(file, ptout,
                                                         std::locale());
           },
           py::arg("file"))
      .def("update",
           [](ComPWA::PartList &partL, ComPWA::ParameterList &pars) {
             UpdateParticleList(partL, pars);
           },
           py::arg("parameters"));

  py::class_<ComPWA::FitParameter, ComPWA::Parameter,
             std::shared_ptr<ComPWA::FitParameter>>(m, "FitParameter")
      .def(py::init<>())
      .def(py::init<std::string, const double, const double>())
      .def("__repr__", &ComPWA::FitParameter::to_str);
  m.def("log", [](const ComPWA::FitParameter p) { LOG(info) << p; });

  py::class_<ComPWA::ParameterList>(m, "ParameterList")
      .def(py::init<>())
      .def("num_parameters", &ComPWA::ParameterList::numParameters)
      .def("__repr__", &ComPWA::ParameterList::to_str);
  m.def("log", [](const ComPWA::ParameterList l) { LOG(info) << l; });

  py::class_<HelicityKinematics, ComPWA::Kinematics,
             std::shared_ptr<HelicityKinematics>>(m, "HelicityKinematics")
      .def(py::init<std::shared_ptr<ComPWA::PartList>, std::vector<ComPWA::pid>,
                    std::vector<ComPWA::pid>, std::array<double, 4>>())
      .def(py::init<std::shared_ptr<ComPWA::PartList>, std::vector<ComPWA::pid>,
                    std::vector<ComPWA::pid>>())
      .def("set_phsp_volume", &ComPWA::Physics::HelicityFormalism::
                                  HelicityKinematics::setPhspVolume)
      .def("print_sub_systems",
           [](const ComPWA::Physics::HelicityFormalism::HelicityKinematics
                  &kin) {
             LOG(info) << "Subsystems used by HelicityKinematics:";
             for (auto i : kin.subSystems()) {
               // Have to add " " here (bug in boost 1.59)
               LOG(info) << " " << i;
             }
           });

  py::class_<ComPWA::Tools::RootGenerator, ComPWA::Generator,
             std::shared_ptr<ComPWA::Tools::RootGenerator>>(m, "RootGenerator")
      .def(py::init<std::shared_ptr<ComPWA::PartList>,
                    std::shared_ptr<ComPWA::Kinematics>, int>());

  py::class_<ComPWA::Tools::RootPlot, std::shared_ptr<ComPWA::Tools::RootPlot>>(
      m, "RootPlot")
      .def(py::init<std::shared_ptr<ComPWA::Kinematics>>())
      .def("set_data", (void (ComPWA::Tools::RootPlot::*)(
                           std::shared_ptr<ComPWA::DataReader::Data>)) &
                           ComPWA::Tools::RootPlot::setDataSample)
      .def("set_phsp_sample", (void (ComPWA::Tools::RootPlot::*)(
                                  std::shared_ptr<ComPWA::DataReader::Data>)) &
                                  ComPWA::Tools::RootPlot::setPhspSample)
      .def("set_intensity", &ComPWA::Tools::RootPlot::setIntensity)
      .def("add_component",&ComPWA::Tools::RootPlot::addComponent)
      .def("write", &ComPWA::Tools::RootPlot::write);

  py::class_<ComPWA::Estimator::MinLogLH, ComPWA::IEstimator,
             std::shared_ptr<ComPWA::Estimator::MinLogLH>>(m, "MinLogLH")
      .def(py::init<std::shared_ptr<ComPWA::Kinematics>,
                    std::shared_ptr<ComPWA::AmpIntensity>,
                    std::shared_ptr<ComPWA::DataReader::Data>,
                    std::shared_ptr<ComPWA::DataReader::Data>,
                    std::shared_ptr<ComPWA::DataReader::Data>, unsigned int,
                    unsigned int>())
      .def("enable_function_tree",
           &ComPWA::Estimator::MinLogLH::UseFunctionTree)
      .def("log_function_tree",
           [](ComPWA::Estimator::MinLogLH &min) {
             LOG(info) << min.tree()->head()->print(25);
           })
      .def("print_function_tree", [](ComPWA::Estimator::MinLogLH &min) {
        return min.tree()->head()->print(25);
      });

  py::class_<ComPWA::Optimizer::Minuit2::MinuitIF, ComPWA::Optimizer::Optimizer,
             std::shared_ptr<ComPWA::Optimizer::Minuit2::MinuitIF>>(m,
                                                                    "MinuitIF")
      .def(py::init<std::shared_ptr<ComPWA::IEstimator>,
                    ComPWA::ParameterList &>())
      .def("enable_hesse", &ComPWA::Optimizer::Minuit2::MinuitIF::setUseHesse)
      .def("minimize", &ComPWA::Optimizer::Minuit2::MinuitIF::exec);

  py::class_<ComPWA::Optimizer::Minuit2::MinuitResult, ComPWA::FitResult,
             std::shared_ptr<ComPWA::Optimizer::Minuit2::MinuitResult>>(
      m, "MinuitResult")
      .def("set_fit_fractions",
           &ComPWA::Optimizer::Minuit2::MinuitResult::setFitFractions)
      .def("fit_fractions",
           &ComPWA::Optimizer::Minuit2::MinuitResult::fitFractions)
      .def("print", &ComPWA::Optimizer::Minuit2::MinuitResult::print,
           py::arg("opt") = "")
      .def("write",
           [](const ComPWA::Optimizer::Minuit2::MinuitResult r,
              std::string file) {
             std::ofstream ofs(file);
             boost::archive::xml_oarchive oa(ofs);
             oa << BOOST_SERIALIZATION_NVP(r);
           },
           py::arg("file"));
}
