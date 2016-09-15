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
using namespace DataReader;

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

    //create dummy final state event to initialized the kinematics class
    unsigned int dataSize = 10000;

    std::shared_ptr<Data> data(new RootReader());
    std::shared_ptr<Data> phsp(new RootReader());
    std::shared_ptr<Generator> gen(new Physics::DPKinematics::RootGenerator());

    RunManager run(dataSize, amp, gen);
    run.setGenerator(gen);
    run.setData(data);
    run.setPhspSample(phsp);
    run.generatePhsp(dataSize * 100);
    run.generate(dataSize);
    std::cout << "Data size: " << data->getNEvents() << std::endl;
    data->writeData("data.root", "events");
    phsp->writeData("phspdata.root", "events");

    //===== Plot amplitude
    /*  TH2D hist("plot", "", 200, 0.0, 10.0, 200, 0.0, 10.0);
     TH2D plot_onlyweight("plot", "", 100, 0.0, 10.0, 100, 0.0, 10.0);

     Event event;

     unsigned int num_events(200000);
     progressBar bar(num_events);
     for (unsigned int i = 0; i < num_events; i++) {
     gen->generate(event);
     double evWeight = event.getWeight();
     ComPWA::dataPoint point(event);
     ParameterList ampPar = amp->intensity(point);
     const ComPWA::Particle& p0 = event.getParticle(0);
     const ComPWA::Particle& p1 = event.getParticle(1);
     const ComPWA::Particle& p2 = event.getParticle(2);

     hist.Fill(p1.invariantMass(p2), p0.invariantMass(p1),
     evWeight * ampPar.GetDoubleParameter(0)->GetValue());
     plot_onlyweight.Fill(p1.invariantMass(p2), p0.invariantMass(p1),
     evWeight);
     }

     hist.Scale(1.0 * num_events / hist.Integral());
     plot_onlyweight.Scale(1.0 * num_events / plot_onlyweight.Integral());

     TCanvas c;*/

    /*// get default styles from xml file
     NeatPlotting::XMLStyleConfigParser xml_style_config_parser(
     "Tools/RootNeatPlotting/default-style.xml");
     NeatPlotting::DefaultStyleSingleton::Instance().readStyleConfigParser(
     xml_style_config_parser);

     // create an empty plot bundle
     NeatPlotting::PlotBundle plot_bundle;

     //create a drawable root object and style pair (default style is being used)
     NeatPlotting::DrawableDataObjectStylePair<TH1*> hist_style_pair(&hist);
     // change the draw option
     hist_style_pair.draw_style.draw_option = "COLZ";

     // and add it to the plot bundle
     // the plot bundle will not make a copy of the histogram
     // so make sure the lifetime of hist is long enough
     plot_bundle.addHistogram(hist_style_pair);

     // add labels to the x and y axis
     plot_bundle.plot_axis.x_axis_title = "M^{2}_#pi_0#pi_0";
     plot_bundle.plot_axis.y_axis_title = "M^{2}_#pi_0#gamma";

     plot_bundle.drawOnCurrentPad();

     c.SetLogz(1);

     hist.Draw("colz");
     c.SaveAs("plot.pdf");
     plot_onlyweight.Draw("colz");
     c.SaveAs("plot_onlyweight.pdf");*/

    //  TFile output(outFile.c_str(),"update");
    //  output.SetCompressionLevel(1); //try level 2 also
    //  output.Close(); //clean output file
    //  plotData plot(std::string("dalitz"),outFile, data);  plot.plot(300);
    //  plotData plotPhsp(std::string("dalitz_phsp"),outFile, data);  plotPhsp.plot(300);
  }
  return 0;
}
