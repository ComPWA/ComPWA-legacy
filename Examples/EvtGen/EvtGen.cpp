

// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Simple Dalitz plot analysis with ComPWA
///

#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/Logging.hpp"
#include "Core/Properties.hpp"
#include "Physics/EvtGen/ComPWAEvtGenIF.hpp"
#include "Physics/EvtGen/DalitzKinematics.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IncoherentIntensity.hpp"
#include "Physics/ParticleList.hpp"
#include "Tools/DalitzPlot.hpp"
#include "Tools/FitFractions.hpp"
#include "Tools/Generate.hpp"
#include "Tools/ParameterTools.hpp"
#include "Tools/RootGenerator.hpp"

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

// Root header files go here
#include "TFile.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLorentzVector.h"

using namespace ComPWA;
using ComPWA::Optimizer::Minuit2::MinuitResult;
using ComPWA::Physics::IncoherentIntensity;
using ComPWA::Physics::EvtGenIF::DalitzKinematics;
using ComPWA::Physics::EvtGenIF::EvtGenIF;
using ComPWA::Physics::HelicityFormalism::HelicityKinematics;

// Enable serialization of MinuitResult. For some reason has to be outside
// any namespaces.
BOOST_CLASS_EXPORT(ComPWA::Optimizer::Minuit2::MinuitResult)

// We define an intensity model using a raw string literal. Currently, this is
// just a toy model without any physical meaning.
// (comments within the string are ignored!). This is convenient since we
// do not have to configure the build system to copy input files somewhere.
// In practise you may want to use a normal XML input file instead.
std::string amplitudeModel = R"####(
<Intensity Class='Incoherent' Name="jpsiGammaPiPi_inc">
  <Intensity Class='Coherent' Name="jpsiGammaPiPi">
    <Amplitude Class="SequentialPartialAmplitude" Name="f2(1270)">
      <Parameter Class='Double' Type="Magnitude"  Name="Magnitude_f2">
        <Value>1.0</Value>
        <Min>-1.0</Min>
        <Max>2.0</Max>
        <Fix>false</Fix>
      </Parameter>
      <Parameter Class='Double' Type="Phase" Name="Phase_f2">
        <Value>0.0</Value>
        <Min>-100</Min>
        <Max>100</Max>
        <Fix>false</Fix>
      </Parameter>
      <PartialAmplitude Class="HelicityDecay" Name="f2ToPiPi">
        <DecayParticle Name="f2(1270)" Helicity="0"/>
        <RecoilSystem FinalState="0" />
        <DecayProducts>
          <Particle Name="pi0" FinalState="1"  Helicity="0"/>
          <Particle Name="pi0" FinalState="2"  Helicity="0"/>
        </DecayProducts>
      </PartialAmplitude>
    </Amplitude>
    <Amplitude Class="SequentialPartialAmplitude" Name="myAmp">
      <Parameter Class='Double' Type="Magnitude"  Name="Magnitude_my">
        <Value>1.0</Value>
        <Min>-1.0</Min>
        <Max>2.0</Max>
        <Fix>true</Fix>
      </Parameter>
      <Parameter Class='Double' Type="Phase" Name="Phase_my`">
        <Value>0.0</Value>
        <Min>-100</Min>
        <Max>100</Max>
        <Fix>true</Fix>
      </Parameter>
      <PartialAmplitude Class="HelicityDecay" Name="MyResToPiPi">
        <DecayParticle Name="myRes" Helicity="0"/>
        <RecoilSystem FinalState="0" />
        <DecayProducts>
          <Particle Name="pi0" FinalState="1"  Helicity="0"/>
          <Particle Name="pi0" FinalState="2"  Helicity="0"/>
        </DecayProducts>
      </PartialAmplitude>
    </Amplitude>
  </Intensity>
</Intensity>
)####";

std::string myParticles = R"####(
<ParticleList>
  <Particle Name="f2(1270)">
    <Pid>225</Pid>
    <Parameter Class='Double' Type="Mass" Name="Mass_f2(1270)">
      <Value>1.2755</Value>
      <Error>8.0E-04</Error>
      <Min>0.1</Min>
      <Max>2.0</Max>
      <Fix>false</Fix>
    </Parameter>
    <QuantumNumber Class="Spin" Type="Spin" Value="2"/>
    <QuantumNumber Class="Int" Type="Charge" Value="0"/>
    <QuantumNumber Class="Int" Type="Parity" Value="+1"/>
    <QuantumNumber Class="Int" Type="Cparity" Value="+1"/>
    <DecayInfo Type="relativisticBreitWigner">
      <FormFactor Type="0" />
      <Parameter Class='Double' Type="Width" Name="Width_f2(1270)">
        <Value>0.1867</Value>
      </Parameter>
      <Parameter Class='Double' Type="MesonRadius" Name="Radius_rho">
        <Value>2.5</Value>
        <Fix>true</Fix>
      </Parameter>
    </DecayInfo>
  </Particle>
  <Particle Name="myRes">
    <Pid>999999</Pid>
    <Parameter Class='Double' Type="Mass" Name="Mass_myRes">
      <Value>2.0</Value>
      <Error>8.0E-04</Error>
    </Parameter>
    <QuantumNumber Class="Spin" Type="Spin" Value="1"/>
    <QuantumNumber Class="Int" Type="Charge" Value="0"/>
    <QuantumNumber Class="Int" Type="Parity" Value="+1"/>
    <QuantumNumber Class="Int" Type="Cparity" Value="+1"/>
    <DecayInfo Type="relativisticBreitWigner">
      <FormFactor Type="0" />
      <Parameter Class='Double' Type="Width" Name="Width_myRes">
        <Value>0.1</Value>
        <Min>0.1</Min>
        <Max>1.0</Max>
        <Fix>false</Fix>
      </Parameter>
      <Parameter Class='Double' Type="MesonRadius" Name="Radius_myRes">
        <Value>2.5</Value>
        <Fix>true</Fix>
      </Parameter>
    </DecayInfo>
  </Particle>
</ParticleList>
)####";

unsigned int nResos = 2;
double resoParameter[2][3] = {{1.2755, 0.186, 2}, {2.0, 0.1, 1}};

///
/// Simple Dalitz plot fit of the channel J/psi -> gamma pi0 pi0
///
/// The basic procedure is the following:
/// 1) Create Kinematics object
/// 2) Generate a large phase space sample
/// 3) Create intensity from pre-defined model
/// 4) Generate a data sample given intensity and kinematics
/// 5) Fit the model to the data
/// 6) Plot data sample and intensity
///
int main(int argc, char **argv) {

  // initialize logging
  Logging log("EvtFit-log.txt", "debug");

  // List with all particle information needed
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, defaultParticleList);
  ReadParticles(partL, myParticles);

  //---------------------------------------------------
  // 1) Create Kinematics object
  //---------------------------------------------------
  std::vector<pid> initialState = {443};
  std::vector<pid> finalState = {22, 111, 111};
  auto kin =
      std::make_shared<DalitzKinematics>(partL, initialState, finalState);

  auto kinB =
      std::make_shared<HelicityKinematics>(partL, initialState, finalState);

  //---------------------------------------------------
  // 2) Generate a large phase space sample
  //---------------------------------------------------
  auto gen = std::make_shared<ComPWA::Tools::RootGenerator>(partL, kin, 173);
  std::shared_ptr<ComPWA::Data::Data> phspSample(
      ComPWA::Tools::generatePhsp(10000000, gen));

  LOG(INFO) << "Phsp generated!";

  //---------------------------------------------------
  // 3) Create intensity from pre-defined model
  //---------------------------------------------------
  // Read in model property_tree
  std::stringstream modelStream;
  modelStream << amplitudeModel;
  boost::property_tree::ptree modelTree;
  boost::property_tree::xml_parser::read_xml(modelStream, modelTree);

  LOG(INFO) << "Model data loaded!";

  double mA = 0., mB = 0.135, mC = 0.135, bigM = 3.0969; // TODO!
  // Construct intensity class from model string
  auto intens = std::make_shared<EvtGenIF>(partL, mA, mB, mC, bigM);
  LOG(INFO) << "EvtGenIF created!";
  intens->addResonances(modelTree.get_child("Intensity"), kin, partL);
  // std::vector<int> recoilS; std::vector<int> finalA; std::vector<int> finalB;
  // recoilS.push_back(3); finalA.push_back(1); finalB.push_back(2);
  // ComPWA::SubSystem subsys(recoilS, finalA, finalB);
  // for(unsigned int i=0; i<nResos; i++){
  //  intens->addResonance("TestBW", resoParameter[i][0], resoParameter[i][1],
  //  resoParameter[i][2], subsys);
  //}
  LOG(INFO) << "Resonances added to EvtGenDalitz";

  // Pass phsp sample to intensity for normalization.
  // Convert to dataPoints first.
  auto phspPoints =
      std::make_shared<std::vector<DataPoint>>(phspSample->dataPoints(kin));
  intens->setPhspSample(phspPoints, phspPoints);

  //---------------------------------------------------
  // 3b) Create intensity from comparisson model
  //---------------------------------------------------

  // Construct intensity class from model string
  auto intensB = std::make_shared<IncoherentIntensity>(
      partL, kinB, modelTree.get_child("Intensity"));

  // Pass phsp sample to intensity for normalization.
  // Convert to dataPoints first.
  auto phspPointsB =
      std::make_shared<std::vector<DataPoint>>(phspSample->dataPoints(kinB));
  intensB->setPhspSample(phspPointsB, phspPointsB);

  //---------------------------------------------------
  // 4) Generate a data sample given intensity and kinematics
  //---------------------------------------------------
  std::shared_ptr<ComPWA::Data::Data> sample = ComPWA::Tools::generate(
      1000000, kin, gen, intens, phspSample, phspSample);
  LOG(INFO) << "EvtGenDalitz events generated";

  std::shared_ptr<ComPWA::Data::Data> sampleB = ComPWA::Tools::generate(
      1000000, kinB, gen, intensB, phspSample, phspSample);
  LOG(INFO) << "EvtGenDalitz events generated";

  //---------------------------------------------------
  // 5) Fit the model to the data and print the result
  //---------------------------------------------------
  ParameterList fitPar;
  intens->parameters(fitPar);
  // Set start error of 0.05 for parameters, run Minos?
  setErrorOnParameterList(fitPar, 0.05, false);
  fitPar.doubleParameter(1)->fixParameter(true);
  fitPar.doubleParameter(3)->fixParameter(true);

  auto esti = std::make_shared<Estimator::MinLogLH>(
      kin, intens, sample, phspSample, phspSample, 0, 0);

  esti->UseFunctionTree(false);
  // esti->tree()->parameter();
  // LOG(debug) << esti->tree()->head()->print(25);

  auto minuitif = new Optimizer::Minuit2::MinuitIF(esti, fitPar);
  minuitif->setUseHesse(false);

  // STARTING MINIMIZATION
  auto result = std::dynamic_pointer_cast<MinuitResult>(minuitif->exec(fitPar));

  // Calculate fit fractions
  // std::vector<std::pair<std::string, std::string>> fitComponents;
  // fitComponents.push_back(
  //     std::pair<std::string, std::string>("myAmp", "jpsiGammaPiPi"));
  // fitComponents.push_back(
  //     std::pair<std::string, std::string>("f2(1270)", "jpsiGammaPiPi"));

  //  ParameterList fitFracs = Tools::CalculateFitFractions(
  //    kin, intens->component("jpsiGammaPiPi"), phspPoints, fitComponents);
  // A proper calculation of the fit fraction uncertainty requires
  // the uncertainty of the fit parameters propagated. We do this
  // using a numerical approach. Using the covariance matrix
  // 100 independend sets of fit parameters are generated and the fit fractions
  // are recalculated. In the end we take the RMS.
  // Tools::CalcFractionError(fitPar, result->covarianceMatrix(), fitFracs, kin,
  //                        intens->component("jpsiGammaPiPi"), phspPoints, 100,
  //                        fitComponents);

  // result->setFitFractions(fitFracs);
  // result->print();

  //---------------------------------------------------
  // 5.1) Save the fit result
  //---------------------------------------------------
  // std::ofstream ofs("DalitzFit-fitResult.xml");
  // boost::archive::xml_oarchive oa(ofs);
  // oa << BOOST_SERIALIZATION_NVP(result);

  // UpdateParticleList(partL, fitPar);
  //  boost::property_tree::ptree ptout;
  //  ptout.add_child("ParticleList", SaveParticles(partL));
  // ptout.add_child("IncoherentIntensity", intens->save());
  // boost::property_tree::xml_parser::write_xml("DalitzFit-Model.xml", ptout,
  //                                           std::locale());
  //---------------------------------------------------
  // 6) Plot data sample and intensity
  //---------------------------------------------------
  // ComPWA::Tools::DalitzPlot pl(kin, "DalitzFit", 100);
  //  pl.setData(sample);
  // pl.setPhspData(phspSample);
  // pl.setFitAmp(intens, "", kBlue - 4);
  // pl.plot();
  LOG(INFO) << "Done";

  std::shared_ptr<ComPWA::Data::Data> sampleFit = ComPWA::Tools::generate(
      1000000, kin, gen, intens, phspSample, phspSample);
  LOG(INFO) << "EvtGenDalitz fitresult events generated";

  unsigned int nBins = 100;
  unsigned int maxEvents = sample->numEvents();
  double masssq12, masssq13, masssq23;
  TH2D *bw12 = new TH2D("bw12", "inv. mass-sq of particles 1&2", nBins, 0., 10.,
                        nBins, 0., 10.);
  bw12->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
  bw12->GetYaxis()->SetTitle("m_{13}^{2} / GeV^{2}/c^{2}");
  bw12->GetXaxis()->CenterTitle();
  bw12->GetYaxis()->CenterTitle();
  TH2D *bw13 = new TH2D("bw13", "inv. mass-sq of particles 1&3", nBins, 0., 10.,
                        nBins, 0., 10.);
  bw13->GetXaxis()->SetTitle("m_{13}^{2} / GeV^{2}/c^{2}");
  bw13->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
  bw13->GetXaxis()->CenterTitle();
  bw13->GetYaxis()->CenterTitle();
  TH2D *bw23 = new TH2D("bw23", "inv. mass-sq of particles 2&3", nBins, 0., 10.,
                        nBins, 0., 10.);
  bw23->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}/c^{2}");
  bw23->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
  bw23->GetXaxis()->CenterTitle();
  bw23->GetYaxis()->CenterTitle();

  unsigned int maxEventsfit = sampleFit->numEvents();
  TH2D *bw12fit = new TH2D("bw12fit", "inv. mass-sq of particles 1&2", nBins,
                           0., 10., nBins, 0., 10.);
  bw12fit->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
  bw12fit->GetYaxis()->SetTitle("m_{13}^{2} / GeV^{2}/c^{2}");
  bw12fit->GetXaxis()->CenterTitle();
  bw12fit->GetYaxis()->CenterTitle();
  TH2D *bw13fit = new TH2D("bw13fit", "inv. mass-sq of particles 1&3", nBins,
                           0., 10., nBins, 0., 10.);
  bw13fit->GetXaxis()->SetTitle("m_{13}^{2} / GeV^{2}/c^{2}");
  bw13fit->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
  bw13fit->GetXaxis()->CenterTitle();
  bw13fit->GetYaxis()->CenterTitle();
  TH2D *bw23fit = new TH2D("bw23fit", "inv. mass-sq of particles 2&3", nBins,
                           0., 10., nBins, 0., 10.);
  bw23fit->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}/c^{2}");
  bw23fit->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
  bw23fit->GetXaxis()->CenterTitle();
  bw23fit->GetYaxis()->CenterTitle();

  unsigned int maxEventsB = sampleB->numEvents();
  // double masssq12B, masssq13B, masssq23B;
  TH2D *bw12B = new TH2D("bw12B", "inv. mass-sq of particles 1&2 Heli", nBins,
                         0., 10., nBins, 0., 10.);
  bw12B->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
  bw12B->GetYaxis()->SetTitle("m_{13}^{2} / GeV^{2}/c^{2}");
  bw12B->GetXaxis()->CenterTitle();
  bw12B->GetYaxis()->CenterTitle();
  TH2D *bw13B = new TH2D("bw13B", "inv. mass-sq of particles 1&3 Heli", nBins,
                         0., 10., nBins, 0., 10.);
  bw13B->GetXaxis()->SetTitle("m_{13}^{2} / GeV^{2}/c^{2}");
  bw13B->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
  bw13B->GetXaxis()->CenterTitle();
  bw13B->GetYaxis()->CenterTitle();
  TH2D *bw23B = new TH2D("bw23B", "inv. mass-sq of particles 2&3 Heli", nBins,
                         0., 10., nBins, 0., 10.);
  bw23B->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}/c^{2}");
  bw23B->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
  bw23B->GetXaxis()->CenterTitle();
  bw23B->GetYaxis()->CenterTitle();

  unsigned int maxEventsPHSP = phspSample->numEvents();
  // double masssq12PHSP, masssq13PHSP, masssq23PHSP;
  TH2D *bw12PHSP = new TH2D("bw12GEN", "inv. mass-sq of particles 1&2 GEN",
                            nBins, 0., 10., nBins, 0., 10.);
  bw12PHSP->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
  bw12PHSP->GetYaxis()->SetTitle("m_{13}^{2} / GeV^{2}/c^{2}");
  bw12PHSP->GetXaxis()->CenterTitle();
  bw12PHSP->GetYaxis()->CenterTitle();
  TH2D *bw13PHSP = new TH2D("bw13GEN", "inv. mass-sq of particles 1&3 GEN",
                            nBins, 0., 10., nBins, 0., 10.);
  bw13PHSP->GetXaxis()->SetTitle("m_{13}^{2} / GeV^{2}/c^{2}");
  bw13PHSP->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
  bw13PHSP->GetXaxis()->CenterTitle();
  bw13PHSP->GetYaxis()->CenterTitle();
  TH2D *bw23PHSP = new TH2D("bw23GEN", "inv. mass-sq of particles 2&3 GEN",
                            nBins, 0., 10., nBins, 0., 10.);
  bw23PHSP->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}/c^{2}");
  bw23PHSP->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}/c^{2}");
  bw23PHSP->GetXaxis()->CenterTitle();
  bw23PHSP->GetYaxis()->CenterTitle();

  for (unsigned int i = 0; i < maxEventsfit; i++) {
    Event event = sampleFit->event(i);

    // myReader.getEvent(-1, a, b, masssq);
    // if(!myReader.getEvent(i, event)) continue; TODO: try exception
    // if(!event.numParticles() == 3) continue;
    // if(!event) continue;
    // cout << "Event: \t" << i << "\t NParticles: \t" << event.getNParticles()
    // << endl;
    const Particle &a(event.particle(0));
    const Particle &b(event.particle(1));
    const Particle &c(event.particle(2));
    masssq12 = pow(a.e() + b.e(), 2) - pow(a.px() + b.px(), 2) -
               pow(a.py() + b.py(), 2) - pow(a.pz() + b.pz(), 2);
    masssq13 = pow(a.e() + c.e(), 2) - pow(a.px() + c.px(), 2) -
               pow(a.py() + c.py(), 2) - pow(a.pz() + c.pz(), 2);
    masssq23 = pow(b.e() + c.e(), 2) - pow(b.px() + c.px(), 2) -
               pow(b.py() + c.py(), 2) - pow(b.pz() + c.pz(), 2);

    bw12fit->Fill(masssq12, masssq13);
    bw13fit->Fill(masssq13, masssq12);
    bw23fit->Fill(masssq23, masssq12);

    // m12->Fill(masssq12);
    // m23->Fill(masssq23);
    // m13->Fill(masssq13);
  }

  for (unsigned int i = 0; i < maxEventsB; i++) {
    Event event = sampleB->event(i);

    // myReader.getEvent(-1, a, b, masssq);
    // if(!myReader.getEvent(i, event)) continue; TODO: try exception
    // if(!event.numParticles() == 3) continue;
    // if(!event) continue;
    // cout << "Event: \t" << i << "\t NParticles: \t" << event.getNParticles()
    // << endl;
    const Particle &a(event.particle(0));
    const Particle &b(event.particle(1));
    const Particle &c(event.particle(2));
    masssq12 = pow(a.e() + b.e(), 2) - pow(a.px() + b.px(), 2) -
               pow(a.py() + b.py(), 2) - pow(a.pz() + b.pz(), 2);
    masssq13 = pow(a.e() + c.e(), 2) - pow(a.px() + c.px(), 2) -
               pow(a.py() + c.py(), 2) - pow(a.pz() + c.pz(), 2);
    masssq23 = pow(b.e() + c.e(), 2) - pow(b.px() + c.px(), 2) -
               pow(b.py() + c.py(), 2) - pow(b.pz() + c.pz(), 2);

    bw12B->Fill(masssq12, masssq13);
    bw13B->Fill(masssq13, masssq12);
    bw23B->Fill(masssq23, masssq12);

    // m12->Fill(masssq12);
    // m23->Fill(masssq23);
    // m13->Fill(masssq13);
  }

  for (unsigned int i = 0; i < maxEvents; i++) {
    Event event = sample->event(i);

    // myReader.getEvent(-1, a, b, masssq);
    // if(!myReader.getEvent(i, event)) continue; TODO: try exception
    if (event.numParticles() != 3)
      continue;
    // if(!event) continue;
    // cout << "Event: \t" << i << "\t NParticles: \t" << event.getNParticles()
    // << endl;
    const Particle &a(event.particle(0));
    const Particle &b(event.particle(1));
    const Particle &c(event.particle(2));
    masssq12 = pow(a.e() + b.e(), 2) - pow(a.px() + b.px(), 2) -
               pow(a.py() + b.py(), 2) - pow(a.pz() + b.pz(), 2);
    masssq13 = pow(a.e() + c.e(), 2) - pow(a.px() + c.px(), 2) -
               pow(a.py() + c.py(), 2) - pow(a.pz() + c.pz(), 2);
    masssq23 = pow(b.e() + c.e(), 2) - pow(b.px() + c.px(), 2) -
               pow(b.py() + c.py(), 2) - pow(b.pz() + c.pz(), 2);

    bw12->Fill(masssq12, masssq13);
    bw13->Fill(masssq13, masssq12);
    bw23->Fill(masssq23, masssq12);

    // m12->Fill(masssq12);
    // m23->Fill(masssq23);
    // m13->Fill(masssq13);
  }

  for (unsigned int i = 0; i < maxEventsPHSP; i++) {
    Event event = phspSample->event(i);

    // myReader.getEvent(-1, a, b, masssq);
    // if(!myReader.getEvent(i, event)) continue; TODO: try exception
    if (event.numParticles() != 3)
      continue;
    // if(!event) continue;
    // cout << "Event: \t" << i << "\t NParticles: \t" << event.getNParticles()
    // << endl;
    const Particle &a(event.particle(0));
    const Particle &b(event.particle(1));
    const Particle &c(event.particle(2));
    masssq12 = pow(a.e() + b.e(), 2) - pow(a.px() + b.px(), 2) -
               pow(a.py() + b.py(), 2) - pow(a.pz() + b.pz(), 2);
    masssq13 = pow(a.e() + c.e(), 2) - pow(a.px() + c.px(), 2) -
               pow(a.py() + c.py(), 2) - pow(a.pz() + c.pz(), 2);
    masssq23 = pow(b.e() + c.e(), 2) - pow(b.px() + c.px(), 2) -
               pow(b.py() + c.py(), 2) - pow(b.pz() + c.pz(), 2);

    bw12PHSP->Fill(masssq12, masssq13);
    bw13PHSP->Fill(masssq13, masssq12);
    bw23PHSP->Fill(masssq23, masssq12);
  }

  TH2D *ratio12 = (TH2D *)bw12->Clone("ratio12");
  TH2D *ratio13 = (TH2D *)bw13->Clone("ratio13");
  TH2D *ratio23 = (TH2D *)bw23->Clone("ratio23");
  ratio12->Divide(bw12B);
  ratio13->Divide(bw13B);
  ratio23->Divide(bw23B);

  TFile output("EvtGenTest.root", "RECREATE", "ROOT_Tree");
  bw12->Write("bw12");
  bw13->Write("bw13");
  bw23->Write("bw23");
  bw12B->Write("bw12Heli");
  bw13B->Write("bw13Heli");
  bw23B->Write("bw23Heli");
  bw12PHSP->Write("phsp12");
  bw13PHSP->Write("phsp13");
  bw23PHSP->Write("phsp23");
  ratio12->Write("12Ratio");
  ratio13->Write("13Ratio");
  ratio23->Write("23Ratio");
  bw12fit->Write("bw12fit");
  bw13fit->Write("bw13fit");
  bw23fit->Write("bw23fit");
  output.Write();
  output.Close();

  std::cout << "Plotting finished " << std::endl;

  return 0;
}
