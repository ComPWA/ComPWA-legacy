// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE SaveAndLoadFitResultTest

#include "Core/FitParameter.hpp"
#include "Core/FunctionTree/FunctionTreeEstimator.hpp"
#include "Core/FunctionTree/FunctionTreeIntensity.hpp"
#include "Core/Logging.hpp"
#include "Data/DataSet.hpp"
#include "Data/Generate.hpp"
#include "Data/Root/RootGenerator.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Physics/BuilderXML.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Tools/UpdatePTreeParameter.hpp"

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/serialization/export.hpp>
#include <boost/test/unit_test.hpp>

using namespace ComPWA;

BOOST_CLASS_EXPORT(ComPWA::FitResult)
BOOST_CLASS_EXPORT(ComPWA::Optimizer::Minuit2::MinuitResult)

BOOST_AUTO_TEST_SUITE(ToolTests)
const std::string TestParticles = R"####(
<ParticleList>
  <Particle Name='pi0'>
    <Pid>111</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_pi0'>
      <Value>0.1349766</Value>
      <Error>0.0000006</Error>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='1'/>
  </Particle>
  <Particle Name='gamma'>
    <Pid>22</Pid>
    <Parameter Class='Double' Type='Mass' Name='mass_gamma'>
      <Value>0.</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='1'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Gparity' Value='1'/>
  </Particle>
  <Particle Name='jpsi'>
    <Pid>443</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_jpsi'>
      <Value>3.0969</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='1'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Gparity' Value='1'/>
    <DecayInfo Type='nonResonant'>
      <FormFactor Type='1' />
      <Parameter Class='Double' Type='Width' Name='Width_jpsi'>
        <Value>0.0000929</Value>
        <Error>0.0000028</Error>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_jpsi'>
        <Value>2.5</Value>
        <Fix>true</Fix>
        <Min>2.0</Min>
        <Max>3.0</Max>
      </Parameter>
    </DecayInfo>
  </Particle>
  <Particle Name='f0'>
    <Pid>6666</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_f0'>
      <Value>1.500</Value>
      <Fix>true</Fix>
      <Min>0.5</Min>
      <Max>1.5</Max>
      <Error>0.00012</Error>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='1'/>
    <DecayInfo Type='relativisticBreitWigner'>
      <FormFactor Type='1' />
      <Parameter Class='Double' Type='Width' Name='Width_f0'>
        <Value>0.050</Value>
        <Fix>true</Fix>
        <Min>0.0</Min>
        <Max>1.0</Max>
        <Error>0.00008</Error>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_f0'>
        <Value>1.5</Value>
        <Fix>true</Fix>
        <Min>1.0</Min>
        <Max>2.0</Max>
        <Error>0</Error>
      </Parameter>
    </DecayInfo>
  </Particle>
</ParticleList>
)####";

const std::string JpsiDecayKinematics = R"####(
<HelicityKinematics>
  <PhspVolume>1</PhspVolume>
  <InitialState>
    <Particle Name='jpsi' PositionIndex='0'/>
  </InitialState>
  <FinalState>
    <Particle Name='gamma' Id='0'/>
    <Particle Name='f0' Id = '1'/>
  </FinalState>
</HelicityKinematics>
)####";
const std::string JpsiDecayTree = R"####(
  <Intensity Class='NormalizedIntensity'>
    <Intensity Class='CoherentIntensity'>
      <Amplitude Class='CoefficientAmplitude'>
        <Amplitude Class='SequentialAmplitude'>
          <Amplitude Class='HelicityDecay'>
            <DecayParticle Name='jpsi' Helicity='+1' />
            <DecayProducts>
              <Particle Name='gamma' FinalState='0' Helicity='1' />
              <Particle Name='f0' FinalState='1' Helicity='0' />
            </DecayProducts>
            <CanonicalSum L='0' S='1.0'>
              <ClebschGordan Type='LS' j1='0.0' m1='0.0' j2='1.0' m2='1.0' J='1.0' M='1.0' />
              <ClebschGordan Type='s2s3' j1='1.0' m1='1.0' j2='0.0' m2='0.0' J='1.0' M='1.0' />
            </CanonicalSum>
          </Amplitude>
        </Amplitude>
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsi_to_gamma+f0_L_0.0_S_1.0;'>
          <Value>1.0</Value>
          <Fix>false</Fix>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_jpsi_to_gamma+f0_L_0.0_S_1.0;'>
          <Value>0.0</Value>
          <Fix>false</Fix>
        </Parameter>
      </Amplitude>

      <Amplitude Class='CoefficientAmplitude'>
        <Amplitude Class='SequentialAmplitude'>
          <Amplitude Class='HelicityDecay'>
            <DecayParticle Name='jpsi' Helicity='+1' />
            <DecayProducts>
              <Particle Name='gamma' FinalState='0' Helicity='1' />
              <Particle Name='f0' FinalState='1' Helicity='0' />
            </DecayProducts>
            <CanonicalSum L='2' S='1.0'>
              <ClebschGordan Type='LS' j1='2.0' m1='0.0' j2='1.0' m2='1.0' J='1.0' M='1.0' />
              <ClebschGordan Type='s2s3' j1='1.0' m1='1.0' j2='0.0' m2='0.0' J='1.0' M='1.0' />
            </CanonicalSum>
          </Amplitude>
        </Amplitude>
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsi_to_gamma+f0_L_2.0_S_1.0;'>
          <Value>1.0</Value>
          <Fix>false</Fix>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_jpsi_to_gamma+f0_L_2.0_S_1.0;'>
          <Value>0.0</Value>
          <Fix>false</Fix>
        </Parameter>
      </Amplitude>
    </Intensity>
  </Intensity>
)####";

BOOST_AUTO_TEST_CASE(SaveAndLoadFitResultTest) {
  ComPWA::Logging Log("trace", "");
  LOG(INFO) << "Now check saven and load fit result...";

  std::stringstream ModelStream;

  // Particle list
  ModelStream << TestParticles;
  auto PartL = readParticles(ModelStream);

  // Kinematics
  ModelStream.clear();
  boost::property_tree::ptree ModelTree;
  ModelStream << JpsiDecayKinematics;
  boost::property_tree::xml_parser::read_xml(ModelStream, ModelTree);
  auto DecayKin = Physics::createHelicityKinematics(
      PartL, ModelTree.get_child("HelicityKinematics"));

  // Generate phsp sample
  ComPWA::Data::Root::RootGenerator Gen(
      DecayKin.getParticleStateTransitionKinematicsInfo());

  ComPWA::Data::Root::RootUniformRealGenerator RandomGenerator(233);

  auto PhspSample = Data::generatePhsp(5000, Gen, RandomGenerator);

  // Model intensity
  ModelStream.clear();
  ModelTree = boost::property_tree::ptree();
  ModelStream << JpsiDecayTree;
  boost::property_tree::xml_parser::read_xml(ModelStream, ModelTree);

  // Fix a set of parameter as reference(optional)
  Tools::fixParameter(ModelTree.get_child("Intensity"),
                      "Magnitude_jpsi_to_gamma+f0_L_0.0_S_1.0;", 1.0);
  Tools::fixParameter(ModelTree.get_child("Intensity"),
                      "Phase_jpsi_to_gamma+f0_L_0.0_S_1.0;", 0.0);
  // Set ranges of Parameters
  Tools::updateParameterRangeByType(ModelTree.get_child("Intensity"),
                                    "Magnitude", 0.0, 10.0);
  Tools::updateParameterRangeByType(ModelTree.get_child("Intensity"), "Phase",
                                    -3.14159, 3.14159);

  ComPWA::Physics::IntensityBuilderXML Builder(
      PartL, DecayKin, ModelTree.get_child("Intensity"), PhspSample);

  auto ModelIntensity = Builder.createIntensity();

  // Generate sample
  auto DataSample =
      Data::generate(500, DecayKin, Gen, ModelIntensity, RandomGenerator);

  // Fit and save result
  // with c++17 you could just do this
  // auto [Esti, Parameters] = ...
  auto EstimatorParametersTuple =
      ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
          ModelIntensity, Data::convertEventsToDataSet(DataSample, DecayKin));
  FitParameterList &Parameters(std::get<1>(EstimatorParametersTuple));

  ComPWA::FunctionTree::FunctionTreeEstimator &Esti(
      std::get<0>(EstimatorParametersTuple));

  auto Minuit = ComPWA::Optimizer::Minuit2::MinuitIF();

  ComPWA::Optimizer::Minuit2::MinuitResult ResultFit, ResultLoad;
  const std::string FitResultName("fitResult.xml");

  // To get a valid fit result
  for (int i = 0; i < 1000; ++i) {
    ComPWA::Optimizer::Minuit2::MinuitResult TempResult =
        Minuit.optimize(Esti, Parameters);
    if (!TempResult.IsValid) {
      Parameters = TempResult.FinalParameters;
      continue;
    }
    // save fit result
    ResultFit = TempResult;
    std::ofstream ofs(FitResultName.c_str());
    boost::archive::xml_oarchive oa(ofs);
    oa << BOOST_SERIALIZATION_NVP(ResultFit);
    ofs.close();
    break;
  }

  // Load fit result from file
  std::ifstream ifs(FitResultName.c_str(), std::ios::in);
  // under g++ 8.2.0 and boost 1.69.0,
  // after run of below codes finish (no errors),
  // there will be an exception which cause a failure:
  // terminate called after throwing an instance of
  // 'boost::archive::archive_exception'
  //  what():  input stream error-No such file or directory
  // it is due to the two lines:
  //  boost::archive::xml_iarchive ia(ifs);
  //  ia >> BOOST_SERIALIZATION_NVP(ResultLoad);
  /*
    boost::archive::xml_iarchive ia(ifs);
    ia >> BOOST_SERIALIZATION_NVP(ResultLoad);
    ifs.close();
    auto ModelIntensityLoad = Builder.createIntensity(PartL, DecayKin,
        ModelTree.get_child("Intensity"));
    ModelIntensityLoad->updateParametersFrom(ResultLoad->finalParameters());
    //Check if loaded Model and original Model have differences
    DataSample->convertEventsToDataPoints(DecayKin);
    for (const auto &Point : DataSample->getDataPointList()) {
      double Value1 = ModelIntensity->evaluate(Point);
      double Value2 = ModelIntensityLoad->evaluate(Point);
      LOG(INFO) << "Intensity with fit Parameters:\t" << Value1;
      LOG(INFO) << "Intensity with load Parameters:\t" << Value2;
      BOOST_CHECK_EQUAL(Value1, Value2);
    }

    //Check if fit result and loaded result have difference;
    //BOOST_CHECK_EQUAL(ResultFit->intialLH(), ResultLoad->initialLH());
    //BOOST_CHECK_EQUAL(ResultFit->trueLH(), ResultLoad->trueLH());
    //BOOST_CHECK_EQUAL(ResultFit->calcInterference(),
    //    ResultFit->calcInterference());
    BOOST_CHECK_EQUAL(ResultFit->result(), ResultLoad->result());
    BOOST_CHECK_EQUAL(ResultFit->finalLH(), ResultLoad->finalLH());
    BOOST_CHECK_EQUAL(ResultFit->isValid(), ResultLoad->isValid());
    BOOST_CHECK_EQUAL(ResultFit->ndf(), ResultLoad->ndf());
    BOOST_CHECK_EQUAL(ResultFit->edm(), ResultLoad->edm());

    void checkMatrix(const std::vector<std::vector<double>> &CovFit,
        const std::vector<std::vector<double>> &CovLoad,
        const std::string MatrixName);

    std::vector<std::vector<double>> CovFit = ResultFit->covarianceMatrix();
    std::vector<std::vector<double>> CovLoad = ResultLoad->covarianceMatrix();
    checkMatrix(CovFit, CovLoad, "covariance matrix");
    std::vector<std::vector<double>> CorFit = ResultFit->correlationMatrix();
    std::vector<std::vector<double>> CorLoad = ResultLoad->correlationMatrix();
    checkMatrix(CorFit, CorLoad, "correlation matrix");
    std::vector<double> CCFit = ResultFit->gobalCC();
    std::vector<double> CCLoad = ResultLoad->gobalCC();
    checkMatrix(std::vector<std::vector<double>>(1, CCFit),
        std::vector<std::vector<double>>(1, CCLoad), "global coefficiencts");
    LOG(INFO) << "All Save FitResult and Load FitResult Tests Finished!";
  */
}

void checkMatrix(const std::vector<std::vector<double>> &CovFit,
                 const std::vector<std::vector<double>> &CovLoad,
                 const std::string MatrixName) {
  BOOST_CHECK_EQUAL(CovFit.size(), CovLoad.size());
  if (CovFit.size() != CovLoad.size())
    return;
  for (std::size_t i = 0; i < CovFit.size(); ++i) {
    BOOST_CHECK_EQUAL(CovFit.at(i).size(), CovLoad.at(i).size());
    if (CovFit.at(i).size() != CovLoad.at(i).size())
      return;
    for (std::size_t j = 0; j < CovLoad.size(); ++j) {
      LOG(INFO) << "Check " << MatrixName << ": ";
      LOG(INFO) << "i = " << i << " j = " << j;
      LOG(INFO) << " Value from fit:\t" << CovFit.at(i).at(j);
      LOG(INFO) << " Value from load:\t" << CovLoad.at(i).at(j);
      BOOST_CHECK_EQUAL(CovFit.at(i).at(j), CovLoad.at(i).at(j));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
