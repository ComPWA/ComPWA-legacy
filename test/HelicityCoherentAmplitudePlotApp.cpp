//! Application to generate J/Psi -> g pi pi.
/*!
 * @file HelicityCoherentAmplitudeGenApp.cpp
 * This application uses the HelicityAmplitude module and a
 * phase-space generator to generate a file with J/Psi -> gamma pi0 pi0 events.
 * Also uses the neat plotting tool to create a dalitz plot of the data.
 */

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"

#include "DataReader/RootGenerator/RootGenerator.hpp"
// Physics Interface header files go here
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/RunManager.hpp"
#include "Core/Efficiency.hpp"
#include "DataReader/RootReader/RootReader.hpp"
#include "Core/ProgressBar.hpp"

#include "Physics/DecayTree/DecayConfiguration.hpp"
#include "Physics/DecayTree/DecayXMLConfigReader.hpp"
#include "Physics/DecayTree/DecayTreeFactory.hpp"
#include "Physics/HelicityAmplitude/TopologyAmplitudeFactory.hpp"
#include "Physics/HelicityAmplitude/HelicityKinematics.hpp"
#include "Physics/HelicityAmplitude/CoherentAmplitude.hpp"

/*#include "Tools/RootNeatPlotting/HelperFunctions.h"
 #include "Tools/RootNeatPlotting/plotting/PlotBundle.h"
 #include "Tools/RootNeatPlotting/plotting/Booky.h"
 #include "Tools/RootNeatPlotting/style/DataObjectStyle.h"
 #include "Tools/RootNeatPlotting/style/DefaultStyleSingleton.h"
 #include "Tools/RootNeatPlotting/style/xml-parser/XMLStyleConfigParser.h"*/
//#include "PWA/PlotData.hpp"
using namespace ComPWA;
using DataReader::Data;
using DataReader::RootReader::RootReader;
using DataReader::RootGenerator::RootGenerator;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv) {
  ComPWA::Physics::HelicityFormalism::HelicityKinematics* kin =
      dynamic_cast<ComPWA::Physics::HelicityFormalism::HelicityKinematics*>(ComPWA::Physics::HelicityFormalism::HelicityKinematics::createInstance());

  //load resonances
  std::string input_config_file("Physics/HelicityAmplitude/JPSI_ypipi.xml");
  if (argc > 1) {
    input_config_file = argv[1];
  }

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
            topology_amplitudes));

    std::shared_ptr<Generator> gen(new Physics::DPKinematics::RootGenerator());
    gen->generate(dummy_event);

    ComPWA::dataPoint dummy_data_point;
    kinematics->translateEventToDataPoint(dummy_event, dummy_data_point);

    TH1D ct_pi_pi("cms_costheta_pi0_pi0_system", "", 100, 0, 3.2);
    TH1D ct_pi_pi_hel("hel_costheta_pi0_pi0_system", "", 100, 0, 3.2);
    TH1D phi_pi_pi("cms_phi_pi0_pi0_system", "", 100, -3.2, 3.2);
    TH1D phi_pi_pi_hel("hel_phi_pi0_pi0_system", "", 100, -3.2, 3.2);

    // ok now just randomize theta and phis

    double genMaxVal(0.0);
    double AMPpdf(0.0);

    gRandom = new TRandom3(123456);

    for (unsigned int i = 0; i < 500000; ++i) {
      dummy_data_point.setVal(3, gRandom->Uniform(0, TMath::Pi()));
      dummy_data_point.setVal(4, gRandom->Uniform(-TMath::Pi(), TMath::Pi()));

      dummy_data_point.setVal(8, gRandom->Uniform(0, TMath::Pi()));
      dummy_data_point.setVal(9, gRandom->Uniform(-TMath::Pi(), TMath::Pi()));

      ParameterList list = amp->intensity(dummy_data_point);    //unfortunatly not thread safe
      AMPpdf = *list.GetDoubleParameter(0);
      if (AMPpdf > genMaxVal)
        genMaxVal = 1.2*AMPpdf;
    }

    unsigned int acceptedEvents(0);
    while (acceptedEvents < 100000) {
      dummy_data_point.setVal(3, gRandom->Uniform(0, TMath::Pi()));
      dummy_data_point.setVal(4, gRandom->Uniform(-TMath::Pi(), TMath::Pi()));

      dummy_data_point.setVal(8, gRandom->Uniform(0, TMath::Pi()));
      dummy_data_point.setVal(9, gRandom->Uniform(-TMath::Pi(), TMath::Pi()));

      double ampRnd = gRandom->Uniform(0.0, 1.0) * genMaxVal;
      ParameterList list;
      list = amp->intensity(dummy_data_point);    //unfortunatly not thread safe
      AMPpdf = *list.GetDoubleParameter(0);
      if (genMaxVal < AMPpdf) {
        std::stringstream ss;
        ss << "RunManager::generate: error in HitMiss procedure. "
            << "Maximum value of random number generation smaller then amplitude maximum! "
            << genMaxVal << " < " << AMPpdf;
        throw std::runtime_error(ss.str());
      }
      if (ampRnd > AMPpdf)
        continue;

      ++acceptedEvents;

      ct_pi_pi.Fill(dummy_data_point.getVal(3));
      ct_pi_pi_hel.Fill(dummy_data_point.getVal(8));
      phi_pi_pi.Fill(dummy_data_point.getVal(4));
      phi_pi_pi_hel.Fill(dummy_data_point.getVal(9));
    }

    TCanvas c;
    c.Divide(2, 4);
    c.cd(2);
    ct_pi_pi.Draw();
    ct_pi_pi.GetYaxis()->SetRangeUser(0.0, ct_pi_pi.GetMaximum());
    c.cd(4);
    ct_pi_pi_hel.Draw();
    ct_pi_pi_hel.GetYaxis()->SetRangeUser(0.0, ct_pi_pi_hel.GetMaximum());

    c.cd(6);
    phi_pi_pi.Draw();
    phi_pi_pi.GetYaxis()->SetRangeUser(0.0, phi_pi_pi.GetMaximum());
    c.cd(8);
    phi_pi_pi_hel.Draw();
    phi_pi_pi_hel.GetYaxis()->SetRangeUser(0.0, phi_pi_pi_hel.GetMaximum());
    c.SaveAs("blub.pdf");

  }
  return 0;
}
