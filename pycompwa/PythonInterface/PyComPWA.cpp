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
#include "Estimator/FunctionTreeEstimator.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IncoherentIntensity.hpp"
#include "Physics/IntensityBuilderXML.hpp"
#include "Physics/ParticleList.hpp"
#include "Physics/ParticleStateTransitionKinematicsInfo.hpp"

#include "Tools/EvtGenGenerator.hpp"
#include "Tools/FitFractions.hpp"
#include "Tools/Generate.hpp"
#include "Tools/ParameterTools.hpp"
#include "Tools/Plotting/RootPlotData.hpp"
#include "Tools/RootGenerator.hpp"
#include "Tools/UpdatePTreeParameter.hpp"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(ComPWA::PartList);
PYBIND11_MAKE_OPAQUE(std::vector<ComPWA::Particle>);
PYBIND11_MAKE_OPAQUE(std::vector<ComPWA::Event>);
PYBIND11_MAKE_OPAQUE(std::vector<ComPWA::DataPoint>);
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<ComPWA::FitParameter>>);

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

  py::class_<ComPWA::Parameter, std::shared_ptr<ComPWA::Parameter>>(
      m, "Parameter");

  py::class_<ComPWA::FitParameter, ComPWA::Parameter,
             std::shared_ptr<ComPWA::FitParameter>>(m, "FitParameter")
      .def("__repr__", &ComPWA::FitParameter::to_str)
      .def_property("is_fixed", &ComPWA::FitParameter::isFixed,
                    &ComPWA::FitParameter::fixParameter)
      .def_property("value", &ComPWA::FitParameter::value,
                    &ComPWA::FitParameter::setValue)
      .def_property_readonly("name", &ComPWA::FitParameter::name)
      .def_property_readonly("error", &ComPWA::FitParameter::error);
  m.def("log", [](const ComPWA::FitParameter p) { LOG(INFO) << p; });

  py::class_<std::vector<std::shared_ptr<ComPWA::FitParameter>>>(
      m, "FitParameterVector")
      .def("index",
           [](std::vector<std::shared_ptr<ComPWA::FitParameter>> &v,
              const std::string &name) {
             auto result = std::find_if(
                 v.begin(), v.end(),
                 [&name](std::shared_ptr<ComPWA::FitParameter> p) -> bool {
                   return p->name() == name;
                 });
             if (result != v.end()) {
               return result - v.begin();
             } else {
               throw py::value_error("Parameter name not found!");
             }
           })
      .def("__repr__",
           [](const std::vector<std::shared_ptr<ComPWA::FitParameter>> &v) {
             std::stringstream ss;
             for (unsigned int i = 0; i < v.size(); ++i) {
               ss << "[" << i << "] Name: \"" << v[i]->name()
                  << "\"  Value: " << v[i]->value();
               if (v[i]->isFixed())
                 ss << "  fixed";
               ss << "\n";
             }
             return ss.str();
           })
      .def("__len__",
           [](const std::vector<std::shared_ptr<ComPWA::FitParameter>> &v) {
             return v.size();
           })
      .def("__iter__",
           [](std::vector<std::shared_ptr<ComPWA::FitParameter>> &v) {
             return py::make_iterator(v.begin(), v.end());
           },
           py::keep_alive<0, 1>())
      .def("__getitem__",
           [](std::vector<std::shared_ptr<ComPWA::FitParameter>> &v,
              unsigned int index) { return v.at(index); });

  py::class_<ComPWA::ParameterList>(m, "ParameterList")
      .def(py::init<>())
      .def("__repr__", &ComPWA::ParameterList::to_str)
      .def("get_fit_parameters",
           [](const ComPWA::ParameterList &list) {
             return list.doubleParameters();
           })
      .def("set_parameter_error",
           [](ComPWA::ParameterList &list, double err, bool asymError) {
             setErrorOnParameterList(list, err, asymError);
           },
           "Set error on all fit parameters. @use_asymmetric_errors "
           "triggers the use of MINOS (asymmetic errors).",
           py::arg("error"), py::arg("use_asymmetric_errors"));
  m.def("log", [](const ComPWA::ParameterList l) { LOG(INFO) << l; },
        "Print ParameterList to logging system.");

  // ------- Parameters in ptree
  m.def("update_parameter_range_by_type",
      ComPWA::Tools::updateParameterRangeByType,
      "Update parameters' range of a ptree by parameter type, e.g., Magnitude.",
      py::arg("tree"), py::arg("parameter_type"), py::arg("min"),
      py::arg("max"));
  m.def("update_parameter_range_by_name",
      ComPWA::Tools::updateParameterRangeByName,
      "Update parameters' range of a ptree by parameter name.",
      py::arg("tree"), py::arg("parameter_name"), py::arg("min"),
      py::arg("max"));
  m.def("update_parameter_value",
      ComPWA::Tools::updateParameterValue,
      "Update parameters' value of a ptree by parameter name.",
      py::arg("tree"), py::arg("parameter_name"), py::arg("value"));
  m.def("fix_parameter",
      ComPWA::Tools::fixParameter,
      "Fix parameters current value (to value) of a ptree by parameter name.",
      py::arg("tree"), py::arg("parameter_name"), py::arg("value") = -999);
  m.def("release_parameter",
      ComPWA::Tools::releaseParameter,
      "Release parameters' value (to new value) of a ptree by parameter name.",
      py::arg("tree"), py::arg("parameter_name"), py::arg("value") = -999);
  m.def("update_parameter",
      (void (*)(boost::property_tree::ptree &, const std::string &,
          const std::string &, double, bool, double, double,
          bool, bool, bool) ) &ComPWA::Tools::updateParameter,
      "Update parameters' value, range, fix status, of a ptree.",
      py::arg("tree"), py::arg("key_type"), py::arg("key_value"),
      py::arg("value"), py::arg("fix"), py::arg("min"), py::arg("max"),
      py::arg("update_value"), py::arg("update_fix"), py::arg("update_range"));
  m.def("update_parameter",
      (void (*)(boost::property_tree::ptree &,
          const std::vector<std::shared_ptr<ComPWA::FitParameter>> &) )
          &ComPWA::Tools::updateParameter,
      "Update parameters according input FitParameters.",
      py::arg("tree"), py::arg("fit_parameters"));

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

  py::class_<ComPWA::DataPoint>(m, "DataPoint")
      .def(py::init<>())
      .def("__repr__",
           [](ComPWA::DataPoint &p) {
             std::stringstream ss;
             ss << p;
             return ss.str();
           })
      .def_readonly("kinematic_variable_list",
                    &ComPWA::DataPoint::KinematicVariableList)
      .def_readonly("weight", &ComPWA::DataPoint::Weight);

  m.def("log", [](const ComPWA::DataPoint p) { LOG(INFO) << p; });

  py::bind_vector<std::vector<ComPWA::DataPoint>>(m, "DataPointList");

  py::class_<ComPWA::Data::DataSet, std::shared_ptr<ComPWA::Data::DataSet>>(
      m, "DataSet")
      .def(py::init<const std::vector<ComPWA::Event> &>())
      .def(py::init<const std::vector<ComPWA::DataPoint> &>())
      .def("get_kinematic_variable_names",
           &ComPWA::Data::DataSet::getKinematicVariableNames)
      .def("get_event_list", &ComPWA::Data::DataSet::getEventList)
      .def("get_data_point_list", &ComPWA::Data::DataSet::getDataPointList)
      .def("convert_events_to_datapoints",
           &ComPWA::Data::DataSet::convertEventsToDataPoints,
           "Internally convert the events to data points.",
           py::arg("kinematics"))
      .def("convert_events_to_parameterlist",
           &ComPWA::Data::DataSet::convertEventsToParameterList,
           "Internally convert the events to a horizontal data structure.",
           py::arg("kinematics"))
      .def("add_intensity_weights", &ComPWA::Data::DataSet::addIntensityWeights,
           "Add the intensity values as weights to this data sample.",
           py::arg("intensity"), py::arg("kinematics"));

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

  py::class_<ComPWA::Intensity, std::shared_ptr<ComPWA::Intensity>>(m,
                                                                    "Intensity")
      .def("add_unique_parameters_to",
           &ComPWA::Intensity::addUniqueParametersTo);

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
        (std::shared_ptr<ComPWA::Data::DataSet>(*)(
            unsigned int, std::shared_ptr<ComPWA::Kinematics>,
            std::shared_ptr<ComPWA::Generator>,
            std::shared_ptr<ComPWA::Intensity>)) &
            ComPWA::Tools::generate,
        "Generate sample from an Intensity", py::arg("size"), py::arg("kin"),
        py::arg("gen"), py::arg("intens"));

  m.def("generate",
        (std::shared_ptr<ComPWA::Data::DataSet>(*)(
            unsigned int, std::shared_ptr<ComPWA::Kinematics>,
            std::shared_ptr<ComPWA::Generator>,
            std::shared_ptr<ComPWA::Intensity>,
            std::shared_ptr<ComPWA::Data::DataSet>,
            std::shared_ptr<ComPWA::Data::DataSet>)) &
            ComPWA::Tools::generate,
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

  py::class_<ComPWA::Estimator::Estimator,
             std::shared_ptr<ComPWA::Estimator::Estimator>>(m, "Estimator");

  py::class_<ComPWA::Estimator::MinLogLH, ComPWA::Estimator::Estimator,
             std::shared_ptr<ComPWA::Estimator::MinLogLH>>(m, "MinLogLH")
      .def(py::init([](std::shared_ptr<ComPWA::Intensity> intensity,
                       std::shared_ptr<ComPWA::Data::DataSet> datasample,
                       std::shared_ptr<ComPWA::Data::DataSet> phspdatasample) {
             if (!datasample) {
               LOG(ERROR) << "PyCompWA::MinLogLH(): no data sample given!";
               return std::shared_ptr<ComPWA::Estimator::MinLogLH>(nullptr);
             }
             std::vector<ComPWA::DataPoint> phsppoints;
             if (phspdatasample)
               phsppoints = phspdatasample->getDataPointList();
             return std::make_shared<ComPWA::Estimator::MinLogLH>(
                 intensity, datasample->getDataPointList(), phsppoints);
           }),
           "Create a minimum log likelihood estimator using data points.");

  py::class_<ComPWA::Estimator::FunctionTreeEstimator,
             ComPWA::Estimator::Estimator,
             std::shared_ptr<ComPWA::Estimator::FunctionTreeEstimator>>(
      m, "FunctionTreeEstimator")
      .def("print_function_tree",
           [](std::shared_ptr<ComPWA::Estimator::FunctionTreeEstimator> est,
              int level) { LOG(INFO) << est->print(level); });

  m.def("create_unbinned_log_likelihood_function_tree_estimator",
        &ComPWA::Estimator::createMinLogLHFunctionTreeEstimator,
        py::arg("intensity"), py::arg("datapoints"), py::arg("phsppoints"));

  py::class_<ComPWA::Optimizer::Optimizer,
             std::shared_ptr<ComPWA::Optimizer::Optimizer>>(m, "Optimizer");

  py::class_<ComPWA::Optimizer::Minuit2::MinuitIF, ComPWA::Optimizer::Optimizer,
             std::shared_ptr<ComPWA::Optimizer::Minuit2::MinuitIF>>(m,
                                                                    "MinuitIF")
      .def(py::init<std::shared_ptr<ComPWA::Estimator::Estimator>,
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

  m.def("fit_fractions", &ComPWA::Tools::calculateFitFractions,
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
      py::arg("components") = std::vector<std::string>());

  //------- Plotting

  m.def("create_data_array",
        [](std::shared_ptr<ComPWA::Data::DataSet> DataSample) {
          auto KinVarNames = DataSample->getKinematicVariableNames();
          KinVarNames.push_back("weight");

          std::vector<std::vector<double>> DataArray;
          auto const &DataPoints = DataSample->getDataPointList();
          DataArray.reserve(DataPoints.size());
          for (auto const &x : DataPoints) {
            auto row(x.KinematicVariableList);
            row.push_back(x.Weight);
            DataArray.push_back(row);
          }
          return std::make_pair(KinVarNames, DataArray);
        },
        py::return_value_policy::move);

  m.def("create_fitresult_array",
        [](std::shared_ptr<ComPWA::Intensity> Intensity,
           std::shared_ptr<ComPWA::Data::DataSet> DataSample) {
          auto KinVarNames = DataSample->getKinematicVariableNames();
          KinVarNames.push_back("intensity");
          KinVarNames.push_back("weight");

          std::vector<std::vector<double>> DataArray;
          auto const &DataPoints = DataSample->getDataPointList();
          DataArray.reserve(DataPoints.size());
          for (auto const &x : DataPoints) {
            auto row(x.KinematicVariableList);
            row.push_back(Intensity->evaluate(x));
            row.push_back(x.Weight);
            DataArray.push_back(row);
          }
          return std::make_pair(KinVarNames, DataArray);
        },
        py::return_value_policy::move);

  m.def(
      "create_rootplotdata",
      [](const std::string &filename, std::shared_ptr<ComPWA::Kinematics> kin,
         std::shared_ptr<ComPWA::Data::DataSet> DataSample,
         std::shared_ptr<ComPWA::Data::DataSet> PhspSample,
         std::shared_ptr<ComPWA::Intensity> Intensity,
         std::map<std::string, std::shared_ptr<const ComPWA::Intensity>>
             IntensityComponents,
         std::shared_ptr<ComPWA::Data::DataSet> HitAndMissSample,
         const std::string &option) {
        try {
          auto KinematicsInfo =
              (std::dynamic_pointer_cast<
                   ComPWA::Physics::HelicityFormalism::HelicityKinematics>(kin)
                   ->getParticleStateTransitionKinematicsInfo());
          if (DataSample || (PhspSample && Intensity) || HitAndMissSample) {
            ComPWA::Tools::Plotting::RootPlotData plotdata(KinematicsInfo,
                                                           filename, option);
            if (DataSample) {
              plotdata.writeData(*DataSample.get());
            }
            if (PhspSample && Intensity) {
              plotdata.writeIntensityWeightedPhspSample(
                  *PhspSample.get(), Intensity, IntensityComponents);
            }
            if (HitAndMissSample) {
              plotdata.writeHitMissSample(*HitAndMissSample.get());
            }
          }
        } catch (const std::exception &e) {
          LOG(ERROR) << e.what();
        }
      },
      py::arg("filename"), py::arg("kinematics"),
      py::arg("data_sample") = std::shared_ptr<ComPWA::Data::DataSet>(nullptr),
      py::arg("phsp_sample") = std::shared_ptr<ComPWA::Data::DataSet>(nullptr),
      py::arg("intensity") = std::shared_ptr<ComPWA::Intensity>(nullptr),
      py::arg("intensity_components") =
          std::map<std::string, std::shared_ptr<const ComPWA::Intensity>>(),
      py::arg("hit_and_miss_sample") =
          std::shared_ptr<ComPWA::Data::DataSet>(nullptr),
      py::arg("tfile_option") = "RECREATE");
}
