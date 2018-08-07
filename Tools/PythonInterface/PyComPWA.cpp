// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// ComPWA Python Interface using pybind11
///


#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

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

#include "Core/Generator.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Kinematics.hpp"

#include "DataReader/RootReader/RootReader.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

#include "Physics/ParticleList.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IncoherentIntensity.hpp"


#include "Tools/ParameterTools.hpp"
#include "Tools/RootGenerator.hpp"
#include "Tools/Plotting/RootPlotData.hpp"
#include "Tools/FitFractions.hpp"
#include "Tools/Generate.hpp"

namespace py = pybind11;

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

  auto intens = std::make_shared<ComPWA::Physics::IncoherentIntensity>(
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

//! Helper Class to provide access to raw event data
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
/// >>> import pycompwa as pwa
/// >>> ...
///
PYBIND11_MODULE(pycompwa, m) {
  m.doc() = "ComPWA python interface"; // optional module docstring

  // ------- Logging

  py::class_<ComPWA::Logging, std::shared_ptr<ComPWA::Logging>>(m, "Logging")
      .def(py::init<std::string, std::string>(), "Initialize logging system",
           py::arg("output"), py::arg("log_level"));
  m.def("log", [](std::string msg) { LOG(INFO) << msg; },
        "Write string to logging system.");

  // ------- Parameters

  py::class_<ComPWA::Parameter, std::shared_ptr<ComPWA::Parameter>>(
      m, "Parameter");

  py::class_<ComPWA::FitParameter, ComPWA::Parameter,
             std::shared_ptr<ComPWA::FitParameter>>(m, "FitParameter")
      .def(py::init<>())
      .def(py::init<std::string, const double, const double>(), py::arg("name"),
           py::arg("value"), py::arg("error"))
      .def(py::init<std::string, const double, const double, const double>(),
           py::arg("name"), py::arg("value"), py::arg("min"), py::arg("max"))
      .def(py::init<std::string, const double, const double, const double,
                    const double>(),
           py::arg("name"), py::arg("value"), py::arg("min"), py::arg("max"),
           py::arg("error"))
      .def("__repr__", &ComPWA::FitParameter::to_str);
  m.def("log", [](const ComPWA::FitParameter p) { LOG(INFO) << p; });

  py::class_<ComPWA::ParameterList>(m, "ParameterList")
      .def(py::init<>())
      .def("__repr__", &ComPWA::ParameterList::to_str)
      .def("num_parameters", &ComPWA::ParameterList::numParameters)
      .def("set_parameter_error",
           [](ComPWA::ParameterList &list, double err, bool asymError) {
             setErrorOnParameterList(list, err, asymError);
           },
           "Set error on all fit parameters. @use_asymmetric_errors "
           "triggers the use of MINOS (asymmetic errors).",
           py::arg("error"), py::arg("use_asymmetric_errors"));
  m.def("log", [](const ComPWA::ParameterList l) { LOG(INFO) << l; },
        "Print ParameterList to logging system.");

  // ------- Data

  py::class_<ComPWA::Particle>(m, "Particle")
      .def(py::init<std::array<double, 4>, int, int>(), "", py::arg("p4"),
           py::arg("pid"), py::arg("charge") = 0)
      .def("__repr__",
           [](const ComPWA::Particle &p) {
             std::stringstream ss;
             ss << p;
             return ss.str();
           })
      .def("p4", [](const ComPWA::Particle &p) { return p.fourMomentum()(); });

  py::class_<ComPWA::Event, std::shared_ptr<ComPWA::Event>>(m, "Event")
      .def(py::init<>())
//      .def("__repr__",
//           [](const ComPWA::Event &ev) {
//             std::stringstream ss;
//             ss << ev;
//             return ss.str();
//           })
      .def("num_particles", &ComPWA::Event::numParticles)
      .def("add_particle", &ComPWA::Event::addParticle)
      .def("cms_energy", &ComPWA::Event::cmsEnergy)
      .def("weight", &ComPWA::Event::weight)
      .def("set_weight", &ComPWA::Event::setWeight)
      .def("efficiency", &ComPWA::Event::efficiency)
      .def("set_efficiency", &ComPWA::Event::setEfficiency);

  py::class_<ComPWA::DataReader::Data,
             std::shared_ptr<ComPWA::DataReader::Data>>(m, "Data")
      .def(py::init<>())
      .def("size", &ComPWA::DataReader::Data::numEvents)
      .def("num_events", &ComPWA::DataReader::Data::numEvents)
      .def("add_event", &ComPWA::DataReader::Data::add);

  py::class_<ComPWA::DataReader::RootReader, ComPWA::DataReader::Data,
             std::shared_ptr<ComPWA::DataReader::RootReader>>(m, "RootReader")
      .def(py::init<>(), "Empty RootReader object.")
      .def(py::init<std::string, std::string, int>(),
           "Read ROOT tree from file.", py::arg("input_file"),
           py::arg("tree_name"), py::arg("size") = -1)
      .def("write", &ComPWA::DataReader::RootReader::writeData,
           "Save data as ROOT tree to file.", py::arg("file"),
           py::arg("tree_name"));

  py::class_<ComPWA::DataPoint>(m, "DataPoint")
      .def(py::init<>())
      .def("__repr__",
           [](ComPWA::DataPoint &p) {
             std::stringstream ss;
             ss << p;
             return ss.str();
           })
      .def("size", &ComPWA::DataPoint::size)
      .def("value", &ComPWA::DataPoint::value)
      .def("weight", &ComPWA::DataPoint::weight)
      .def("efficiency", &ComPWA::DataPoint::efficiency);
  m.def("log", [](const ComPWA::DataPoint p) { LOG(INFO) << p; });

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

  // ------- Particles

  py::class_<ComPWA::PartList, std::shared_ptr<ComPWA::PartList>>(m, "PartList")
      .def(py::init<>())
      .def("__repr__",
           [](const ComPWA::PartList &p) {
             std::stringstream ss;
             ss << p;
             return ss.str();
           })
      .def("write",
           [](const ComPWA::PartList &list, std::string file) {
             boost::property_tree::ptree ptout;
             ptout.add_child("ParticleList", SaveParticles(list));
             boost::property_tree::xml_parser::write_xml(file, ptout,
                                                         std::locale());
           },
           "Write particle list to XML.", py::arg("file"))
      .def("update",
           [](ComPWA::PartList &partL, ComPWA::ParameterList &pars) {
             UpdateParticleList(partL, pars);
           },
           "Update particles in list using parameter list",
           py::arg("parameters"));

  m.def("read_particles",
        (void (*)(std::shared_ptr<ComPWA::PartList>, std::string)) &
            ComPWA::ReadParticles);

  m.def("default_particles", []() { return defaultParticleList; },
        "Get a list of predefined particles.");

  // ------- Kinematics

  py::class_<ComPWA::Kinematics, std::shared_ptr<ComPWA::Kinematics>>(
      m, "Kinematics")
      .def("convert", &ComPWA::Kinematics::convert,
           "Convert event to DataPoint.")
      .def("phsp_volume", &ComPWA::Kinematics::phspVolume,
           "Convert event to DataPoint.")
      .def("set_phsp_volume", &ComPWA::Kinematics::setPhspVolume);

  py::class_<
      ComPWA::Physics::HelicityFormalism::HelicityKinematics,
      ComPWA::Kinematics,
      std::shared_ptr<ComPWA::Physics::HelicityFormalism::HelicityKinematics>>(
      m, "HelicityKinematics")
      .def(py::init<std::shared_ptr<ComPWA::PartList>, std::vector<ComPWA::pid>,
                    std::vector<ComPWA::pid>, std::array<double, 4>>())
      .def(py::init<std::shared_ptr<ComPWA::PartList>, std::vector<ComPWA::pid>,
                    std::vector<ComPWA::pid>>())
      .def(py::init([](std::shared_ptr<ComPWA::PartList> partL, std::string model_file_path) {
		    boost::property_tree::ptree pt;
		    boost::property_tree::xml_parser::read_xml(model_file_path, pt);
		    return std::make_shared<ComPWA::Physics::HelicityFormalism::HelicityKinematics>(
				    partL, pt.get_child("HelicityKinematics"));
	    }))
      .def("print_sub_systems",
           [](const ComPWA::Physics::HelicityFormalism::HelicityKinematics
                  &kin) {
             LOG(INFO) << "Subsystems used by HelicityKinematics:";
             for (auto i : kin.subSystems()) {
               // Have to add " " here (bug in boost 1.59)
               LOG(INFO) << " " << i;
             }
           });

  // ------- Intensity

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
           "Write intensity model to XML file.", py::arg("file"));

  m.def("incoherent_intensity",
        (std::shared_ptr<ComPWA::AmpIntensity>(*)(
            std::string, std::shared_ptr<ComPWA::PartList>,
            std::shared_ptr<ComPWA::Kinematics>,
            std::shared_ptr<ComPWA::DataReader::Data>,
            std::shared_ptr<ComPWA::DataReader::Data>)) &
            createIntens,
        "Create an incoherent intensity from a model file. A list of "
        "particles, a kinematics object and one (or two) phase space samples "
        "are needed.",
        py::arg("model_file"), py::arg("particle_list"), py::arg("kin"),
        py::arg("phsp_sample"), py::arg("toy_phsp_sample"));

  //------- Generate

  py::class_<ComPWA::Generator, std::shared_ptr<ComPWA::Generator>>(
      m, "Generator");

  py::class_<ComPWA::Tools::RootGenerator, ComPWA::Generator,
             std::shared_ptr<ComPWA::Tools::RootGenerator>>(m, "RootGenerator")
      .def(py::init<std::shared_ptr<ComPWA::PartList>,
                    std::shared_ptr<ComPWA::Kinematics>, int>());

  m.def("generate", (bool (*)(int, std::shared_ptr<ComPWA::Kinematics>,
                              std::shared_ptr<ComPWA::Generator>,
                              std::shared_ptr<ComPWA::AmpIntensity>,
                              std::shared_ptr<ComPWA::DataReader::Data>,
                              std::shared_ptr<ComPWA::DataReader::Data>,
                              std::shared_ptr<ComPWA::DataReader::Data>)) &
                        ComPWA::Tools::generate,
        "Generate sample from AmpIntensity. In case that detector "
        "reconstruction and selection is considered in the phase space sample "
        "a second pure toy sample needs to be passed. If no phase space "
        "samples are passed events are generated on the fly.",
        py::arg("size"), py::arg("kin"), py::arg("gen"), py::arg("intens"),
        py::arg("sample"),
        py::arg("phspSample") = std::shared_ptr<ComPWA::DataReader::Data>(),
        py::arg("toyPhspSample") = std::shared_ptr<ComPWA::DataReader::Data>());

  m.def("generate_phsp", (bool (*)(int, std::shared_ptr<ComPWA::Generator>,
                                   std::shared_ptr<ComPWA::DataReader::Data>)) &
                             ComPWA::Tools::generatePhsp,
        "Generate phase space sample");

  //------- Estimator + Optimizer

  py::class_<ComPWA::IEstimator, std::shared_ptr<ComPWA::IEstimator>>(
      m, "Estimator");

  py::class_<ComPWA::Estimator::MinLogLH, ComPWA::IEstimator,
             std::shared_ptr<ComPWA::Estimator::MinLogLH>>(m, "MinLogLH")
      .def(py::init<std::shared_ptr<ComPWA::Kinematics>,
                    std::shared_ptr<ComPWA::AmpIntensity>,
                    std::shared_ptr<ComPWA::DataReader::Data>,
                    std::shared_ptr<ComPWA::DataReader::Data>,
                    std::shared_ptr<ComPWA::DataReader::Data>, unsigned int,
                    unsigned int>())
      .def("enable_function_tree",
           &ComPWA::Estimator::MinLogLH::UseFunctionTree,
           "Enable FunctionTree (a caching infrastructure) for the calculation "
           "of the likelihood.")
      .def("log_function_tree",
           [](ComPWA::Estimator::MinLogLH &min) {
             LOG(INFO) << min.tree()->head()->print(25);
           },
           "Print FunctionTree to logging system.")
      .def("print_function_tree",
           [](ComPWA::Estimator::MinLogLH &min) {
             return min.tree()->head()->print(25);
           },
           "Return FunctionTree structure as a string.");

  py::class_<ComPWA::Optimizer::Optimizer,
             std::shared_ptr<ComPWA::Optimizer::Optimizer>>(m, "Optimizer");

  py::class_<ComPWA::Optimizer::Minuit2::MinuitIF, ComPWA::Optimizer::Optimizer,
             std::shared_ptr<ComPWA::Optimizer::Minuit2::MinuitIF>>(m,
                                                                    "MinuitIF")
      .def(py::init<std::shared_ptr<ComPWA::IEstimator>,
                    ComPWA::ParameterList &>())
      .def("enable_hesse", &ComPWA::Optimizer::Minuit2::MinuitIF::setUseHesse,
           "Enable the usage of HESSE after MIGRAD has found a minimum.")
      .def("minimize", &ComPWA::Optimizer::Minuit2::MinuitIF::exec,
           "Start minimization.");

  //------- FitResult

  py::class_<ComPWA::FitResult, std::shared_ptr<ComPWA::FitResult>>(m,
                                                                    "FitResult")
      .def("final_parameters", &ComPWA::FitResult::finalParameters);

  py::class_<ComPWA::Optimizer::Minuit2::MinuitResult, ComPWA::FitResult,
             std::shared_ptr<ComPWA::Optimizer::Minuit2::MinuitResult>>(
      m, "MinuitResult")
      .def("set_fit_fractions",
           &ComPWA::Optimizer::Minuit2::MinuitResult::setFitFractions)
      .def("fit_fractions",
           &ComPWA::Optimizer::Minuit2::MinuitResult::fitFractions)
      .def("log", &ComPWA::Optimizer::Minuit2::MinuitResult::print,
           py::arg("opt") = "", "Print fit result to the logging system.")
      .def("write",
           [](const ComPWA::Optimizer::Minuit2::MinuitResult r,
              std::string file) {
             std::ofstream ofs(file);
             boost::archive::xml_oarchive oa(ofs);
             oa << BOOST_SERIALIZATION_NVP(r);
           },
           py::arg("file"));

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
        "intensity. components is a list of tuples of names e.g.[Amplitude, "
        "AmpIntensity] for which fit fractions are calculated.",
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
        auto resultM =
            std::dynamic_pointer_cast<ComPWA::Optimizer::Minuit2::MinuitResult>(
                fitResult);
        auto phspPoints = std::make_shared<std::vector<ComPWA::DataPoint>>(
            phspSample->dataPoints(kin));
        for (auto i : components)
          std::cout << i.first << " " << i.second << std::endl;
        ComPWA::Tools::CalcFractionError(
            fitParameters, resultM->covarianceMatrix(), fitFractions, kin,
            intens, phspPoints, nSets, components);
      },
      "Calculate uncertainties fot a list of fit fractions given a fit result "
      "and its corresponding intensity.",
      py::arg("fit_params"), py::arg("fit_result"), py::arg("fit_fractions"),
      py::arg("intensity"), py::arg("components"), py::arg("kin"),
      py::arg("phspSample"), py::arg("nSets"));

  //------- Plotting

  py::class_<ComPWA::Tools::Plotting::RootPlotData, std::shared_ptr<ComPWA::Tools::Plotting::RootPlotData>>(
      m, "RootPlotData")
      .def(py::init<std::shared_ptr<ComPWA::Kinematics>, std::shared_ptr<ComPWA::AmpIntensity>>())
      .def("use_efficiency_correction", (void (ComPWA::Tools::Plotting::RootPlotData::*)(
                           std::shared_ptr<ComPWA::DataReader::Data>)) &
                           ComPWA::Tools::Plotting::RootPlotData::useEfficiencyCorrection)
      .def("set_data", (void (ComPWA::Tools::Plotting::RootPlotData::*)(
                           std::shared_ptr<ComPWA::DataReader::Data>)) &
                           ComPWA::Tools::Plotting::RootPlotData::setData)
      .def("set_phsp_mc", (void (ComPWA::Tools::Plotting::RootPlotData::*)(
                           std::shared_ptr<ComPWA::DataReader::Data>)) &
                           ComPWA::Tools::Plotting::RootPlotData::setPhspMC)
      .def("set_hit_miss_mc", (void (ComPWA::Tools::Plotting::RootPlotData::*)(
                           std::shared_ptr<ComPWA::DataReader::Data>)) &
                           ComPWA::Tools::Plotting::RootPlotData::setHitMissMC)
      .def("add_component", &ComPWA::Tools::Plotting::RootPlotData::addComponent,
           "Add component for which weights should be calculated. A tuple of "
           "e.g.[Amplitude, "
           "AmpIntensity] is expected.")
      .def("write", &ComPWA::Tools::Plotting::RootPlotData::write);
}

