//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Test-Application for full fit with simple BW-dalitz-model.
/*!
 * @file DalitzFitApp.cpp
 * This tiny application tests a dalitz-fit procedure with a simple resonance
 * model. It uses the simple LH-estimator MinLogLH, it reads data
 * via the root-reader module RootReader and uses the intensity provided by
 * the Breit-Wigner-Sum  physics module AmplitudeSum. The optimization of the
 * parameters is done with the Minuit2 module MinuitIF. As result the
 * optimized parameters are printed to the terminal.
 */

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

// Root header files go here
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH2D.h"
#include "TMath.h"
#include "TParticle.h"
#include "TGenPhaseSpace.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TRandom3.h"

//Core header files go here
#include "Core/DataPoint.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"

// ComPWA header files go here
#include "DataReader/RootReader/RootReader.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Efficiency.hpp"
#include "Core/ProgressBar.hpp"
#include "Core/DataPointStorage.hpp"

#include "Physics/DecayTree/DecayConfiguration.hpp"
#include "Physics/DecayTree/DecayXMLConfigReader.hpp"
#include "Physics/DecayTree/DecayTreeFactory.hpp"
#include "Physics/HelicityAmplitude/TopologyAmplitudeFactory.hpp"
#include "Physics/HelicityAmplitude/HelicityKinematics.hpp"
#include "Physics/HelicityAmplitude/CoherentAmplitude.hpp"

#include <boost/filesystem.hpp>

const Double_t M = 3.096916;    // GeV/c² (J/psi+)
const Double_t Br = 0.000093;    // GeV/c² (width)
const Double_t m1 = 0.;    // GeV/c² (gamma)
const Double_t m2 = 0.139570;    // GeV/c² (pi)
const Double_t m3 = 0.139570;    // GeV/c² (pi)
//const Double_t c = 299792458.; // m/s
const Double_t PI = 3.14159;    // m/s

using namespace ComPWA;
using DataReader::RootReader::RootReader;
using Estimator::MinLogLH::MinLogLH;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv) {
  int seed = 1234;
  string input_config_file("Physics/HelicityAmplitude/JPSI_ypipi.xml");
  string input_data_file_path("3Part-4vecs.root");
  string input_phsp_file_path("3Part-4vecs.root");
  string output_filename("FitResultJPSI.root");
  string output_directory("./");
  string output_file_suffix("");
  if (argc > 1)
    input_config_file = argv[1];
  if (argc > 2)
    input_data_file_path = argv[2];
  if (argc > 3)
    input_phsp_file_path = argv[3];
  if (argc > 4)
    output_filename = argv[4];
  if (argc > 5)
    output_directory = argv[5];
  if (argc > 6)
    seed = atoi(argv[6]);
  if (argc > 7)
    output_file_suffix = argv[7];

  Logging log("log", boost::log::trivial::debug);    //initialize logging

  BOOST_LOG_TRIVIAL(info)<< "  ComPWA Copyright (C) 2013  Mathias Michel ";
  BOOST_LOG_TRIVIAL(info)<< "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt";
  BOOST_LOG_TRIVIAL(info)<< std::endl;

  bool useFctTree = true, resultGen = true;

  BOOST_LOG_TRIVIAL(info)<< "Load Modules";

  ComPWA::Physics::HelicityFormalism::HelicityKinematics* kin =
      dynamic_cast<ComPWA::Physics::HelicityFormalism::HelicityKinematics*>(ComPWA::Physics::HelicityFormalism::HelicityKinematics::createInstance());

  ComPWA::Physics::DecayTree::DecayConfiguration decay_configuration;
  ComPWA::Physics::DecayTree::DecayXMLConfigReader xml_reader(
      decay_configuration);
  xml_reader.readConfig(input_config_file);

  ComPWA::Physics::DecayTree::DecayTreeFactory decay_tree_factory(
      decay_configuration);

  std::vector<ComPWA::Physics::DecayTree::DecayTree> decay_trees =
      decay_tree_factory.createDecayTrees();

  std::cout << "created " << decay_trees.size() << " decay trees from "
      << input_config_file << " config file!" << std::endl;

  if (decay_trees.size() > 0) {
    ComPWA::Physics::HelicityFormalism::TopologyAmplitudeFactory topology_amp_factory;

    Event dummy_event = topology_amp_factory.createDummyEvent(decay_trees[0]);

    std::vector<ComPWA::Physics::DecayTree::DecayNode> leaves =
        decay_trees[0].getLeaves();
    std::vector<ComPWA::Physics::HelicityFormalism::ParticleStateInfo> fs_particles;
    for (auto iter = leaves.begin(); iter != leaves.end(); ++iter) {
      fs_particles.push_back(iter->state_info_);
    }

    ComPWA::Physics::HelicityFormalism::FinalStateParticleCombinatorics fsp_combinatorics;
    fsp_combinatorics.init(fs_particles, dummy_event);

    ComPWA::Physics::HelicityFormalism::HelicityKinematics* kinematics =
        (ComPWA::Physics::HelicityFormalism::HelicityKinematics*) ComPWA::Physics::HelicityFormalism::HelicityKinematics::createInstance();

    std::vector<ComPWA::Physics::HelicityFormalism::TwoBodyDecayTopology> decay_topologies =
        topology_amp_factory.generateDecayTopologies(decay_trees);

    kinematics->setDecayTopologies(decay_topologies);
    kinematics->init(fsp_combinatorics);

    // The helicity amplitude tree factory sorts the decay trees based
    // on their topology, and then creates the list of amplitude trees from the
    // topology grouped decay trees

    std::vector<ComPWA::Physics::HelicityFormalism::TopologyAmplitude> topology_amplitudes =
        topology_amp_factory.generateTopologyAmplitudes(decay_trees);

    std::cout << "created " << topology_amplitudes.size()
        << " topology amplitudes from the decay trees!" << std::endl;

    std::shared_ptr<ComPWA::Physics::HelicityFormalism::CoherentAmplitude> amp(
        new ComPWA::Physics::HelicityFormalism::CoherentAmplitude(
            topology_amplitudes)
    );


    std::shared_ptr<RootReader> myReader(
        new RootReader(input_data_file_path, "events"));
    myReader->resetWeights();    //setting weights to 1
    std::shared_ptr<RootReader> myPHSPReader(
        new RootReader(input_phsp_file_path, "events"));
    myPHSPReader->setEfficiency(shared_ptr<Efficiency>(new UnitEfficiency()));    //setting efficiency to 1

    if (myReader) {
      Event tmp;
      unsigned int num_events(myReader->getNEvents());
      if (num_events > 0) {
        tmp = myReader->getEvent(0);
        DataPointStorage::Instance().layoutDataStorageStructure(1, num_events,
            tmp);
        progressBar bar(num_events);
        for (unsigned int i = 0; i < num_events; ++i) {
          tmp = myReader->getEvent(i);
          DataPointStorage::Instance().addEvent(1, tmp);
          bar.nextEvent();
        }
      }
    }
    if (myPHSPReader) {
      Event tmp;
      unsigned int num_events(myPHSPReader->getNEvents());
      if (num_events > 0) {
        tmp = myPHSPReader->getEvent(0);
        DataPointStorage::Instance().layoutDataStorageStructure(0, num_events,
            tmp);
        progressBar bar(num_events);
        for (unsigned int i = 0; i < num_events; ++i) {
          tmp = myPHSPReader->getEvent(i);
          DataPointStorage::Instance().addEvent(0, tmp);
          bar.nextEvent();
        }
      }
    }

    ParameterList par;
    std::shared_ptr<Optimizer::ControlParameter> esti;
    amp->copyParameterList(par);    //perfect startvalues

    for (auto const& param : par.GetDoubleParameters()) {
      std::cout << param->GetName() << " " << param->IsFixed() << std::endl;
    }
    esti = MinLogLH::createInstance(amp, myReader, myPHSPReader, 0,
        myReader->getNEvents());
    MinLogLH* contrPar = dynamic_cast<MinLogLH*>(&*(esti->Instance()));
    std::shared_ptr<FunctionTree> tree;
    if (useFctTree) {
      contrPar->setUseFunctionTree(1);
      tree = contrPar->getTree();
    }

    std::shared_ptr<Optimizer::Optimizer> opti(
        new Optimizer::Minuit2::MinuitIF(esti, par));

    ParameterList test;

    BOOST_LOG_TRIVIAL(info)<< "LH with optimal parameters: " << esti->controlParameter(par);
    if (useFctTree)
      BOOST_LOG_TRIVIAL(info)<<tree;
    double startInt[par.GetNDouble()], optiInt[par.GetNDouble()];
    for (unsigned int i = 0; i < par.GetNDouble(); i++) {
      std::shared_ptr<DoubleParameter> tmp = par.GetDoubleParameter(i);
      optiInt[i] = tmp->GetValue();
      if (!tmp->IsFixed()) {
        BOOST_LOG_TRIVIAL(debug)<< *tmp;
        tmp->SetValue(tmp->GetValue()*1.2);
        tmp->SetError(tmp->GetValue());
        if (!tmp->GetValue())
        tmp->SetError(1.);
      }
      startInt[i] = tmp->GetValue();
    }
    BOOST_LOG_TRIVIAL(info)<< "LH with following parameters: " << esti->controlParameter(par);
    for (unsigned int i = 0; i < par.GetNDouble(); i++) {
      BOOST_LOG_TRIVIAL(info)<< par.GetDoubleParameter(i)->GetName() << " = " << par.GetDoubleParameter(i)->GetValue();
    }

    BOOST_LOG_TRIVIAL(info)<< "Start Fit";
    std::shared_ptr<FitResult> genResult = opti->exec(par);
    BOOST_LOG_TRIVIAL(info)<< "Final LH = " << genResult->getResult();

    BOOST_LOG_TRIVIAL(info)<< "Optimierte intensitäten: " << esti->controlParameter(par);
    for (unsigned int i = 0; i < par.GetNDouble(); i++) {
      BOOST_LOG_TRIVIAL(info)<< par.GetDoubleParameter(i)->GetName() << " = " << par.GetDoubleParameter(i)->GetValue()
      << "   [ start: " << startInt[i] << " ," << " optimal: " << optiInt[i] << " ]";
    }

    string output_file_path = output_directory + "/fitresult.txt";
    if (output_file_suffix != "")
      output_file_path = output_directory + "/fitresult_" + output_file_suffix
          + ".txt";
    //genResult->writeText(output_file_path);
    output_file_path = output_directory + "/simplefitresult.txt";
    if (output_file_suffix != "")
      output_file_path = output_directory + "/simplefitresult_"
          + output_file_suffix + ".txt";
    genResult->writeSimpleText(output_file_path);

    if (!resultGen)
      return 0;

    BOOST_LOG_TRIVIAL(info)<< "Done";
  }

  return 0;
}
