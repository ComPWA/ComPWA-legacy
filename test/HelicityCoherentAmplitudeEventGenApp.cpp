//! Application to generate J/Psi -> g pi pi.
/*!
 * @file GenDalitzApp.cpp
 * This application uses the Breit-Wigner-Amplitude-Sum module and a
 * phase-space generator to generate a file with J/Psi -> g pi pi events.
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

#include "DataReader/RootGenerator/RootGenerator.hpp"
// Physics Interface header files go here
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/RunManager.hpp"
#include "Core/Efficiency.hpp"
#include "DataReader/RootReader/RootReader.hpp"
#include "Core/ProgressBar.hpp"

#include "Physics/HelicityAmplitude/DecayConfiguration.hpp"
#include "Physics/HelicityAmplitude/DecayXMLConfigReader.hpp"
#include "Physics/HelicityAmplitude/DecayTreeFactory.hpp"
#include "Physics/HelicityAmplitude/TopologyAmplitudeFactory.hpp"
#include "Physics/HelicityAmplitude/HelicityKinematics.hpp"
#include "Physics/HelicityAmplitude/CoherentAmplitude.hpp"

//#include "PWA/PlotData.hpp"

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv) {
  HelicityFormalism::HelicityKinematics* kin =
      dynamic_cast<HelicityFormalism::HelicityKinematics*>(HelicityFormalism::HelicityKinematics::createInstance());

  //load resonances
  std::string input_config_file("Physics/HelicityAmplitude/JPSI_ypipi.xml");

  HelicityFormalism::DecayConfiguration decay_configuration;
  HelicityFormalism::DecayXMLConfigReader xml_reader(decay_configuration);
  xml_reader.readConfig(input_config_file);

  HelicityFormalism::DecayTreeFactory decay_tree_factory(decay_configuration);

  std::vector<HelicityFormalism::DecayTree> decay_trees =
      decay_tree_factory.createDecayTrees();

  std::cout << "created " << decay_trees.size() << " decay trees from "
      << input_config_file << " config file!" << std::endl;

  HelicityFormalism::TopologyAmplitudeFactory topology_amp_factory;

  HelicityFormalism::HelicityKinematics* kinematics =
      (HelicityFormalism::HelicityKinematics*) HelicityFormalism::HelicityKinematics::createInstance();

  std::vector<HelicityFormalism::TwoBodyDecayTopology> decay_topologies =
      topology_amp_factory.generateDecayTopologies(decay_trees);

  kinematics->setDecayTopologies(decay_topologies);

  Event dummy_event = topology_amp_factory.createDummyEvent(decay_trees[0]);
  kinematics->init(dummy_event);

  // The helicity amplitude tree factory sorts the decay trees based
  // on their topology, and then creates the list of amplitude trees from the
  // topology grouped decay trees

  std::vector<HelicityFormalism::TopologyAmplitude> topology_amplitudes =
      topology_amp_factory.generateTopologyAmplitudes(decay_trees);

  std::cout << "created " << topology_amplitudes.size()
      << " topology amplitudes from the decay trees!" << std::endl;

  std::shared_ptr<Generator> gen(new RootGenerator());

  std::shared_ptr<HelicityFormalism::CoherentAmplitude> amp(
      new HelicityFormalism::CoherentAmplitude(topology_amplitudes));
  amp->init(dummy_event);

  //===== Plot amplitude
  /* TH2D plot("plot", "", 100, 0.0, 10.0, 100, 0.0, 10.0);
   TH2D plot_onlyweight("plot", "", 100, 0.0, 10.0, 100, 0.0, 10.0);

   Event event;

   unsigned int num_events(100000);
   progressBar bar(num_events);
   for (unsigned int i = 0; i < num_events; i++) {
   gen->generate(event);
   double evWeight = event.getWeight();
   dataPoint point(event);
   ParameterList ampPar = amp->intensity(point);
   const Particle& p0 = event.getParticle(0);
   const Particle& p1 = event.getParticle(1);
   const Particle& p2 = event.getParticle(2);

   plot.Fill(p1.invariantMass(p2), p0.invariantMass(p1), evWeight*ampPar.GetDoubleParameter(0)->GetValue());
   plot_onlyweight.Fill(p1.invariantMass(p2), p0.invariantMass(p1), evWeight);
   }

   TCanvas c;
   plot.Draw("colz");
   c.SaveAs("plot.pdf");
   plot_onlyweight.Draw("colz");
   c.SaveAs("plot_onlyweight.pdf");*/

  //create dummy final state event to initialized the kinematics class
  if (decay_trees.size() > 0) {
    std::string outFile = "3Part-4vecs.root";
    unsigned int dataSize = 1000;

    Event dummy_event = topology_amp_factory.createDummyEvent(decay_trees[0]);
    kinematics->init(dummy_event);

    std::shared_ptr<Data> data(new RootReader());
    std::shared_ptr<Data> phsp(new RootReader());
    std::shared_ptr<Generator> gen(new RootGenerator());
    std::shared_ptr<HelicityFormalism::CoherentAmplitude> amp(
        new HelicityFormalism::CoherentAmplitude(topology_amplitudes));
    amp->init(dummy_event);

    RunManager run(dataSize, amp, gen);
    run.setGenerator(gen);
    run.setData(data);
    run.setPhspSample(phsp);
    run.generatePhsp(dataSize * 10);
    run.generate(dataSize);
    std::cout << "Data size: " << data->getNEvents() << std::endl;
    data->writeData("data.root", "events");
    phsp->writeData("phspdata.root", "events");

    //  TFile output(outFile.c_str(),"update");
    //  output.SetCompressionLevel(1); //try level 2 also
    //  output.Close(); //clean output file
    //  plotData plot(std::string("dalitz"),outFile, data);  plot.plot(300);
    //  plotData plotPhsp(std::string("dalitz_phsp"),outFile, data);  plotPhsp.plot(300);
  }
  return 0;
}
