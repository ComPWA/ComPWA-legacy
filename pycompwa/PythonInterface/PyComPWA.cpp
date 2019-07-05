// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include "Core/Event.hpp"
#include "Core/Generator.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Particle.hpp"
#include "Data/DataSet.hpp"
#include "Data/RootIO/RootDataIO.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

#include "Core/FunctionTree/FunctionTreeIntensityWrapper.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IntensityBuilderXML.hpp"
#include "Physics/ParticleList.hpp"
#include "Physics/ParticleStateTransitionKinematicsInfo.hpp"

#include "Tools/EvtGenGenerator.hpp"
#include "Tools/FitFractions.hpp"
#include "Tools/Generate.hpp"
#include "Tools/Plotting/RootPlotData.hpp"
#include "Tools/RootGenerator.hpp"
#include "Tools/UpdatePTreeParameter.hpp"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(ComPWA::PartList);
PYBIND11_MAKE_OPAQUE(std::vector<ComPWA::Particle>);
PYBIND11_MAKE_OPAQUE(std::vector<ComPWA::Event>);
PYBIND11_MAKE_OPAQUE(std::vector<ComPWA::DataPoint>);
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(ui, m) {
  m.doc() = "pycompwa module\n"
            "---------------\n";
  // ------- Logging
  // reinitialize the logger with level INFO
  ComPWA::Logging("INFO");
  py::class_<ComPWA::Logging, std::shared_ptr<ComPWA::Logging>>(m, "Logging")
      .def(py::init<std::string, std::string>(), "Initialize logging system",
           py::arg("log_level"), py::arg("filename"))
      .def(py::init<std::string>(), "Initialize logging system",
           py::arg("log_level"));
  m.def("log", [](std::string msg) { LOG(INFO) << msg; },
        "Write string to logging system.");

  // ------- Parameters

  py::class_<ComPWA::FitParameter<double>>(m, "FitParameter")
      .def("__repr__",
           [](const ComPWA::FitParameter<double> &x) {
             std::stringstream ss;
             ss << x;
             return ss.str();
           })
      .def_readwrite("is_fixed", &ComPWA::FitParameter<double>::IsFixed)
      .def_readwrite("value", &ComPWA::FitParameter<double>::Value)
      .def_readwrite("name", &ComPWA::FitParameter<double>::Name)
      .def_readwrite("error", &ComPWA::FitParameter<double>::Error)
      .def_readwrite("bounds", &ComPWA::FitParameter<double>::Bounds);
  m.def("log", [](const ComPWA::FitParameter<double> p) { LOG(INFO) << p; });

  py::class_<ComPWA::FitParameterList>(m, "FitParameterList");

  m.def("log",
        [](const ComPWA::FitParameterList &list) {
          for (auto x : list)
            LOG(INFO) << x;
        },
        "Print FitParameter list to logging system.");

  // ------- Parameters in ptree
  m.def(
      "update_parameter_range_by_type",
      ComPWA::Tools::updateParameterRangeByType,
      "Update parameters' range of a ptree by parameter type, e.g., Magnitude.",
      py::arg("tree"), py::arg("parameter_type"), py::arg("min"),
      py::arg("max"));
  m.def("update_parameter_range_by_name",
        ComPWA::Tools::updateParameterRangeByName,
        "Update parameters' range of a ptree by parameter name.",
        py::arg("tree"), py::arg("parameter_name"), py::arg("min"),
        py::arg("max"));
  m.def("update_parameter_value", ComPWA::Tools::updateParameterValue,
        "Update parameters' value of a ptree by parameter name.",
        py::arg("tree"), py::arg("parameter_name"), py::arg("value"));
  m.def("fix_parameter", ComPWA::Tools::fixParameter,
        "Fix parameters current value (to value) of a ptree by parameter name.",
        py::arg("tree"), py::arg("parameter_name"), py::arg("value") = -999);
  m.def(
      "release_parameter", ComPWA::Tools::releaseParameter,
      "Release parameters' value (to new value) of a ptree by parameter name.",
      py::arg("tree"), py::arg("parameter_name"), py::arg("value") = -999);
  m.def("update_parameter",
        (void (*)(boost::property_tree::ptree &, const std::string &,
                  const std::string &, double, bool, double, double, bool, bool,
                  bool)) &
            ComPWA::Tools::updateParameter,
        "Update parameters' value, range, fix status, of a ptree.",
        py::arg("tree"), py::arg("key_type"), py::arg("key_value"),
        py::arg("value"), py::arg("fix"), py::arg("min"), py::arg("max"),
        py::arg("update_value"), py::arg("update_fix"),
        py::arg("update_range"));
  m.def("update_parameter",
        (void (*)(boost::property_tree::ptree &,
                  const ComPWA::FitParameterList &)) &
            ComPWA::Tools::updateParameter,
        "Update parameters according input FitParameters.", py::arg("tree"),
        py::arg("fit_parameters"));

  // ------- Data

  py::class_<ComPWA::Particle>(m, "Particle")
      .def(py::init<std::array<double, 4>, int>(), "", py::arg("p4"),
           py::arg("pid"))
      .def("__repr__",
           [](const ComPWA::Particle &p) {
             std::stringstream ss;
             ss << p;
             return ss.str();
           })
      .def("p4", [](const ComPWA::Particle &p) { return p.fourMomentum()(); });

  py::bind_vector<std::vector<ComPWA::Particle>>(m, "ParticleList");
  py::class_<ComPWA::Event>(m, "Event")
      .def(py::init<>())
      .def("particle_list",
           [](const ComPWA::Event &ev) { return ev.ParticleList; })
      .def("weight", [](const ComPWA::Event &ev) { return ev.Weight; });

  py::bind_vector<std::vector<ComPWA::Event>>(m, "EventList");

  py::class_<ComPWA::Data::RootDataIO,
             std::shared_ptr<ComPWA::Data::RootDataIO>>(m, "RootDataIO")
      .def(py::init<const std::string &, int>())
      .def(py::init<const std::string &>())
      .def(py::init<>())
      .def("readData", &ComPWA::Data::RootDataIO::readData,
           "Read ROOT tree from file.", py::arg("input_file"))
      .def("writeData", &ComPWA::Data::RootDataIO::writeData,
           "Save data as ROOT tree to file.", py::arg("data"), py::arg("file"));

  m.def("log", [](const ComPWA::DataPoint p) { LOG(INFO) << p; });

  py::class_<ComPWA::Data::DataSet>(m, "DataSet")
      .def_readonly("data", &ComPWA::Data::DataSet::Data)
      .def_readonly("weights", &ComPWA::Data::DataSet::Weights)
      .def_readonly("variable_names", &ComPWA::Data::DataSet::VariableNames);

  m.def("convert_events_to_dataset", &ComPWA::Data::convertEventsToDataSet,
        "Internally convert the events to data points.", py::arg("events"),
        py::arg("kinematics"));
  m.def("add_intensity_weights", &ComPWA::Data::addIntensityWeights,
        "Add the intensity values as weights to this data sample.",
        py::arg("intensity"), py::arg("events"), py::arg("kinematics"));

  // ------- Particles
  py::class_<ComPWA::PartList, std::shared_ptr<ComPWA::PartList>>(m, "PartList")
      .def(py::init<>())
      .def("__repr__", [](const ComPWA::PartList &p) {
        std::stringstream ss;
        ss << p;
        return ss.str();
      });

  m.def("read_particles",
        (void (*)(std::shared_ptr<ComPWA::PartList>, std::string)) &
            ComPWA::ReadParticles);

  m.def("default_particles",
        []() { return ComPWA::Physics::defaultParticleList; },
        "Get a list of predefined particles.");

  // ------- Kinematics

  py::class_<ComPWA::Kinematics, std::shared_ptr<ComPWA::Kinematics>>(
      m, "Kinematics")
      .def("convert", &ComPWA::Kinematics::convert,
           "Convert event to DataPoint.")
      .def("get_kinematic_variable_names",
           &ComPWA::Kinematics::getKinematicVariableNames)
      .def("phsp_volume", &ComPWA::Kinematics::phspVolume,
           "Convert event to DataPoint.");

  py::class_<ComPWA::Physics::ParticleStateTransitionKinematicsInfo>(
      m, "ParticleStateTransitionKinematicsInfo");

  py::class_<
      ComPWA::Physics::HelicityFormalism::HelicityKinematics,
      ComPWA::Kinematics,
      std::shared_ptr<ComPWA::Physics::HelicityFormalism::HelicityKinematics>>(
      m, "HelicityKinematics")
      .def(py::init<std::shared_ptr<ComPWA::PartList>, std::vector<ComPWA::pid>,
                    std::vector<ComPWA::pid>, std::array<double, 4>>())
      .def(py::init<std::shared_ptr<ComPWA::PartList>, std::vector<ComPWA::pid>,
                    std::vector<ComPWA::pid>>())
      .def(py::init<
           const ComPWA::Physics::ParticleStateTransitionKinematicsInfo &,
           double>())
      .def(py::init<
           const ComPWA::Physics::ParticleStateTransitionKinematicsInfo &>())
      .def("create_all_subsystems", &ComPWA::Physics::HelicityFormalism::
                                        HelicityKinematics::createAllSubsystems)
      .def("get_particle_state_transition_kinematics_info",
           &ComPWA::Physics::HelicityFormalism::HelicityKinematics::
               getParticleStateTransitionKinematicsInfo)
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

  py::class_<ComPWA::Intensity, std::shared_ptr<ComPWA::Intensity>>(
      m, "Intensity");

  py::class_<
      ComPWA::FunctionTree::FunctionTreeIntensityWrapper, ComPWA::Intensity,
      std::shared_ptr<ComPWA::FunctionTree::FunctionTreeIntensityWrapper>>(
      m, "FunctionTreeIntensityWrapper")
      .def("evaluate",
           &ComPWA::FunctionTree::FunctionTreeIntensityWrapper::evaluate)
      .def("updateParametersFrom",
           [](ComPWA::FunctionTree::FunctionTreeIntensityWrapper &x,
              ComPWA::FitParameterList pars) {
             std::vector<double> params;
             for (auto x : pars)
               params.push_back(x.Value);
             x.updateParametersFrom(params);
           });

  m.def(
      "create_intensity_and_kinematics",
      [&](const std::string &filename) {
        boost::property_tree::ptree pt;
        boost::property_tree::xml_parser::read_xml(filename, pt);
        ComPWA::Physics::IntensityBuilderXML Builder;
        return Builder.createIntensityAndKinematics(pt);
      },
      "Create an intensity and a helicity kinematics from a xml file. The file "
      "should contain a particle list, and a kinematics and intensity section.",
      py::arg("xml_filename"));

  //------- Generate

  py::class_<ComPWA::Generator, std::shared_ptr<ComPWA::Generator>>(
      m, "Generator");

  py::class_<ComPWA::Tools::RootGenerator, ComPWA::Generator,
             std::shared_ptr<ComPWA::Tools::RootGenerator>>(m, "RootGenerator")
      .def(py::init<
           const ComPWA::Physics::ParticleStateTransitionKinematicsInfo &,
           int>());

  py::class_<ComPWA::Tools::EvtGenGenerator, ComPWA::Generator,
             std::shared_ptr<ComPWA::Tools::EvtGenGenerator>>(m,
                                                              "EvtGenGenerator")
      .def(py::init<
           const ComPWA::Physics::ParticleStateTransitionKinematicsInfo &,
           unsigned int>());

  m.def("generate",
        [](unsigned int n, std::shared_ptr<ComPWA::Kinematics> kin,
           std::shared_ptr<ComPWA::Generator> gen,
           std::shared_ptr<ComPWA::Intensity> intens) {
          return ComPWA::Tools::generate(n, kin, gen, intens);
        },
        "Generate sample from an Intensity", py::arg("size"), py::arg("kin"),
        py::arg("gen"), py::arg("intens"));

  m.def("generate",
        [](unsigned int n, std::shared_ptr<ComPWA::Kinematics> kin,
           std::shared_ptr<ComPWA::Generator> gen,
           std::shared_ptr<ComPWA::Intensity> intens,
           const std::vector<ComPWA::Event> &phspsample) {
          return ComPWA::Tools::generate(n, kin, gen, intens, phspsample);
        },
        "Generate sample from an Intensity, using a given phase space sample.",
        py::arg("size"), py::arg("kin"), py::arg("gen"), py::arg("intens"),
        py::arg("phspSample"));

  m.def("generate",
        [](unsigned int n, std::shared_ptr<ComPWA::Kinematics> kin,
           std::shared_ptr<ComPWA::Generator> gen,
           std::shared_ptr<ComPWA::Intensity> intens,
           const std::vector<ComPWA::Event> &phspsample,
           const std::vector<ComPWA::Event> &toyphspsample) {
          return ComPWA::Tools::generate(n, kin, gen, intens, phspsample,
                                         toyphspsample);
        },
        "Generate sample from an Intensity. In case that detector "
        "reconstruction and selection is considered in the phase space sample "
        "a second pure toy sample needs to be passed.",
        py::arg("size"), py::arg("kin"), py::arg("gen"), py::arg("intens"),
        py::arg("phspSample"), py::arg("toyPhspSample") = nullptr);

  m.def("generate_phsp", &ComPWA::Tools::generatePhsp,
        "Generate phase space sample");

  m.def("generate_importance_sampled_phsp",
        &ComPWA::Tools::generateImportanceSampledPhsp,
        "Generate an Intensity importance weighted phase space sample",
        py::arg("size"), py::arg("kin"), py::arg("gen"), py::arg("intens"));

  //------- Estimator + Optimizer

  py::class_<ComPWA::Estimator::Estimator<double>>(m, "Estimator");

  py::class_<ComPWA::FunctionTree::FunctionTreeEstimatorWrapper,
             ComPWA::Estimator::Estimator<double>>(
      m, "FunctionTreeEstimatorWrapper");

  m.def("create_unbinned_log_likelihood_function_tree_estimator",
        (std::tuple<ComPWA::FunctionTree::FunctionTreeEstimatorWrapper,
                    ComPWA::FitParameterList>(*)(
            std::shared_ptr<ComPWA::FunctionTree::FunctionTreeIntensityWrapper>,
            const ComPWA::Data::DataSet &)) &
            ComPWA::Estimator::createMinLogLHFunctionTreeEstimator,
        py::arg("intensity"), py::arg("datapoints"));

  m.def("create_unbinned_log_likelihood_function_tree_estimator",
        (std::tuple<ComPWA::FunctionTree::FunctionTreeEstimatorWrapper,
                    ComPWA::FitParameterList>(*)(
            std::shared_ptr<ComPWA::FunctionTree::FunctionTreeIntensityWrapper>,
            const ComPWA::Data::DataSet &, const ComPWA::Data::DataSet &)) &
            ComPWA::Estimator::createMinLogLHFunctionTreeEstimator,
        py::arg("intensity"), py::arg("datapoints"), py::arg("phsppoints"));

  py::class_<
      ComPWA::Optimizer::Optimizer<ComPWA::Optimizer::Minuit2::MinuitResult>>(
      m, "Optimizer");

  py::class_<
      ComPWA::Optimizer::Minuit2::MinuitIF,
      ComPWA::Optimizer::Optimizer<ComPWA::Optimizer::Minuit2::MinuitResult>>(
      m, "MinuitIF")
      .def(py::init<>())
      .def("optimize", &ComPWA::Optimizer::Minuit2::MinuitIF::optimize,
           "Start minimization.");

  //------- FitResult

  py::class_<ComPWA::FitResult, std::shared_ptr<ComPWA::FitResult>>(m,
                                                                    "FitResult")
      .def_readonly("final_parameters", &ComPWA::FitResult::FinalParameters)
      .def_readonly("initial_parameters", &ComPWA::FitResult::InitialParameters)
      .def_readonly("initial_estimator_value",
                    &ComPWA::FitResult::InitialEstimatorValue)
      .def_readonly("final_estimator_value",
                    &ComPWA::FitResult::FinalEstimatorValue)
      .def_property_readonly(
          "fit_duration_in_seconds",
          [](const ComPWA::FitResult &x) { return x.FitDuration.count(); })
      .def_readonly("covariance_matrix", &ComPWA::FitResult::CovarianceMatrix);

  py::class_<ComPWA::Optimizer::Minuit2::MinuitResult, ComPWA::FitResult,
             std::shared_ptr<ComPWA::Optimizer::Minuit2::MinuitResult>>(
      m, "MinuitResult")
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

  /*m.def("fit_fractions", &ComPWA::Tools::calculateFitFractions,
        "Calculates the fit fractions for all components of a given coherent "
        "intensity.",
        py::arg("intensity"), py::arg("sample"),
        py::arg("components") = std::vector<std::string>());

  m.def(
      "fit_fractions_with_propagated_errors",
      [](std::shared_ptr<const ComPWA::Physics::CoherentIntensity> CohIntensity,
         std::shared_ptr<ComPWA::Data::DataSet> Sample,
         std::shared_ptr<ComPWA::Optimizer::Minuit2::MinuitResult> Result,
         const std::vector<std::string> &Components) {
        ComPWA::Tools::calculateFitFractionsWithCovarianceErrorPropagation(
            CohIntensity, Sample, Result->covarianceMatrix(), Components);
      },
      "Calculates the fit fractions and errors for all components of a given "
      "coherent intensity.",
      py::arg("intensity"), py::arg("sample"), py::arg("fit_result"),
      py::arg("components") = std::vector<std::string>());*/

  //------- Plotting

  m.def("create_data_array",
        [](ComPWA::Data::DataSet DataSample) {
          auto KinVarNames = DataSample.VariableNames;
          KinVarNames.push_back("weight");

          std::vector<std::vector<double>> DataArray(DataSample.Data);
          DataArray.push_back(DataSample.Weights);
          return std::make_pair(KinVarNames, DataArray);
        },
        py::return_value_policy::move);

  m.def("create_fitresult_array",
        [](std::shared_ptr<ComPWA::Intensity> Intensity,
           ComPWA::Data::DataSet DataSample) {
          auto KinVarNames = DataSample.VariableNames;
          KinVarNames.push_back("intensity");
          KinVarNames.push_back("weight");

          std::vector<std::vector<double>> DataArray(DataSample.Data);
          DataArray.push_back(DataSample.Weights);
          DataArray.push_back(Intensity->evaluate(DataSample.Data));
          return std::make_pair(KinVarNames, DataArray);
        },
        py::return_value_policy::move);

  m.def(
      "create_rootplotdata",
      [](const std::string &filename, std::shared_ptr<ComPWA::Kinematics> kin,
         const ComPWA::Data::DataSet &DataSample,
         const ComPWA::Data::DataSet &PhspSample,
         std::shared_ptr<ComPWA::Intensity> Intensity,
         std::map<std::string, std::shared_ptr<ComPWA::Intensity>>
             IntensityComponents,
         const ComPWA::Data::DataSet &HitAndMissSample,
         const std::string &option) {
        try {
          auto KinematicsInfo =
              (std::dynamic_pointer_cast<
                   ComPWA::Physics::HelicityFormalism::HelicityKinematics>(kin)
                   ->getParticleStateTransitionKinematicsInfo());
          ComPWA::Tools::Plotting::RootPlotData plotdata(KinematicsInfo,
                                                         filename, option);
          plotdata.writeData(DataSample);
          if (Intensity) {
            plotdata.writeIntensityWeightedPhspSample(PhspSample, Intensity,
                                                      IntensityComponents);
          }
          plotdata.writeHitMissSample(HitAndMissSample);
        } catch (const std::exception &e) {
          LOG(ERROR) << e.what();
        }
      },
      py::arg("filename"), py::arg("kinematics"), py::arg("data_sample"),
      py::arg("phsp_sample") = ComPWA::Data::DataSet(),
      py::arg("intensity") = std::shared_ptr<ComPWA::Intensity>(nullptr),
      py::arg("intensity_components") =
          std::map<std::string, std::shared_ptr<ComPWA::Intensity>>(),
      py::arg("hit_and_miss_sample") = ComPWA::Data::DataSet(),
      py::arg("tfile_option") = "RECREATE");
}
