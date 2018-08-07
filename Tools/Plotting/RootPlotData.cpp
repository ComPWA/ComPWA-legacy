// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <stdio.h>
#include <numeric>

#include "Core/DataPoint.hpp"
#include "Core/Logging.hpp"
#include "Core/ProgressBar.hpp"
<<<<<<< HEAD
#include "Core/Properties.hpp"
=======
>>>>>>> 417d84dd2... - added RootPlotData class based on RootPlot class.
#include "Physics/PartialAmplitude.hpp"

#include "Tools/Plotting/RootPlotData.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

namespace ComPWA {
namespace Tools {
namespace Plotting {

RootPlotData::RootPlotData(std::shared_ptr<ComPWA::Kinematics> kin,
    std::shared_ptr<ComPWA::AmpIntensity> intens) :
		Kinematics(kin), Intensity(intens), CorrectForEfficiency(false) {
}

RootPlotData::~RootPlotData() {
}

/// Add sub component of the Intensity. For each event in the phase space
/// sample each component is evaluated and its value is added to the TTree.
void RootPlotData::addComponent(std::string componentName, std::string intensityName, std::string title) {
	std::string ComponentLabel(title);
	if (ComponentLabel == "") ComponentLabel = componentName;

	std::shared_ptr<ComPWA::AmpIntensity> comp;
	try {
		comp = Intensity->component(intensityName)->component(componentName);
	}
	catch (std::exception &ex) {
<<<<<<< HEAD
		LOG(ERROR)<< "RootPlotData::addComponent() | Component " << componentName
=======
		LOG(error)<< "RootPlotData::addComponent() | Component " << componentName
>>>>>>> 417d84dd2... - added RootPlotData class based on RootPlot class.
		<< " of " << componentName <<" not found in AmpIntensity "
		<< Intensity->name() << ".";
		return;
	}
	AmplitudeComponents[ComponentLabel] = comp;
}

<<<<<<< HEAD
void RootPlotData::write(std::string treePrefix, std::string fileName,
                         std::string option) {

  TFile *tf = new TFile(fileName.c_str(), option.c_str());

  auto varNames = Kinematics->variableNames();

  // write out a final state id to name name mapping
  const auto& KinProps(Kinematics->getKinematicsProperties());
  for (unsigned int i = 0; i < KinProps.FinalState.size(); ++i) {
	unsigned int FinalStateID(KinProps.FinalStateEventPositionMapping[i]);
    TObjString s(FindParticle(KinProps.ParticleList, KinProps.FinalState[i]).name().c_str());
    s.Write("" + FinalStateID);
  }

  double EventWeight(1.0);
  double EventEfficiency(1.0);

  double DataIntegral(0.0);
  //===== write data
  if (Data && Data->numEvents() > 0) {
    TTree *tree = new TTree((treePrefix + "_data").c_str(), "DataSample");
    auto DataPointValues = std::vector<double>(varNames.size(), 0.0);
    for (unsigned int i = 0; i < varNames.size(); ++i)
      tree->Branch(varNames.at(i).c_str(), &DataPointValues.at(i),
                   (varNames.at(i) + "/D").c_str());
    tree->Branch("event_weight", &EventWeight, "event_weight/D");
    tree->Branch("event_efficiency", &EventEfficiency, "event_efficiency/D");

    ComPWA::ProgressBar bar(Data->numEvents());
    for (unsigned int i = 0; i < Data->numEvents(); ++i) { // loop over data
      Event event(Data->event(i));

      EventEfficiency = 1.0;
      if (CorrectForEfficiency)
        EventEfficiency = event.efficiency();
      if (EventEfficiency == 0.0) {
        LOG(ERROR) << "DalitzPlot::Fill() | Loop over "
                      "data sample: An event with zero efficiency was found! "
                      "This should not happen! We skip it!";
        continue;
      }
      EventWeight = event.weight();

      /* Construct DataPoint from Event to check if dataPoint is within the
       * phase space boundaries */
      DataPoint point;
      try {
        Kinematics->convert(event, point);
      } catch (BeyondPhsp &ex) { // event outside phase, remove
        continue;
      }

      // Fill branch references with dataPoint
      for (unsigned int j = 0; j < DataPointValues.size(); ++j) {
        DataPointValues[j] = point.value(j);
      }

      DataIntegral += point.weight();
      tree->Fill();
      bar.next();
    }
    tree->Write();
    delete (tree);
  }

  double PhspIntegral(0.0);
  //===== write amplitude weights
  if (WeightedPhspMC && WeightedPhspMC->numEvents() > 0) {
    TTree *tree = new TTree((treePrefix + "_weighted_phsp_MC").c_str(),
                            "WeightedPhspMCSample");
    auto DataPointValues = std::vector<double>(varNames.size(), 0.0);
    for (unsigned int i = 0; i < varNames.size(); ++i)
      tree->Branch(varNames.at(i).c_str(), &DataPointValues.at(i),
                   (varNames.at(i) + "/D").c_str());
    tree->Branch("event_weight", &EventWeight, "event_weight/D");
    tree->Branch("event_efficiency", &EventEfficiency, "event_efficiency/D");

    double IntensityWeight(0.0);
    tree->Branch("intensity", &IntensityWeight, "intensity/D");
    std::vector<double> AmplitudeComponentWeights =
        std::vector<double>(AmplitudeComponents.size(), 0.0);
    unsigned int counter(0);
    for (auto const &amp : AmplitudeComponents) {
      tree->Branch(amp.first.c_str(), &AmplitudeComponentWeights.at(counter),
                   (amp.first + "/D").c_str());
      ++counter;
    }

    ComPWA::ProgressBar bar(WeightedPhspMC->numEvents());
    for (unsigned int i = 0; i < WeightedPhspMC->numEvents();
         i++) { // loop over data
      Event event(WeightedPhspMC->event(i));

      EventEfficiency = 1.0;
      if (CorrectForEfficiency)
        EventEfficiency = event.efficiency();
      if (EventEfficiency == 0.0) {
        LOG(ERROR) << "DalitzPlot::Fill() | Loop over "
                      "data sample: An event with zero efficiency was found! "
                      "This should not happen! We skip it!";
        continue;
      }
      EventWeight = event.weight();
      PhspIntegral += EventWeight;
      /* Construct DataPoint from Event to check if dataPoint is within the
       * phase space boundaries */
      DataPoint point;
      try {
        Kinematics->convert(event, point);
      } catch (BeyondPhsp &ex) { // event outside phase, remove
        continue;
      }

      // Fill branch references with dataPoint
      for (unsigned int j = 0; j < DataPointValues.size(); ++j) {
        DataPointValues[j] = point.value(j);
      }

      IntensityWeight = Intensity->intensity(point);
      // Loop over all components that we want to plot
      counter = 0;
      for (auto const &amp : AmplitudeComponents) {
        AmplitudeComponentWeights[counter] = amp.second->intensity(point);
        ++counter;
      }

      tree->Fill();
      bar.next();
    }
    tree->Write();
    delete (tree);
  }

  //===== write data
  if (HitAndMissMC && HitAndMissMC->numEvents() > 0) {
    TTree *tree =
        new TTree((treePrefix + "_hitmiss_MC").c_str(), "Hit&MissMCSample");
    auto DataPointValues = std::vector<double>(varNames.size(), 0.0);
    for (unsigned int i = 0; i < varNames.size(); ++i)
      tree->Branch(varNames.at(i).c_str(), &DataPointValues.at(i),
                   (varNames.at(i) + "/D").c_str());
    tree->Branch("event_weight", &EventWeight, "event_weight/D");
    tree->Branch("event_efficiency", &EventEfficiency, "event_efficiency/D");

    ComPWA::ProgressBar bar(HitAndMissMC->numEvents());
    for (unsigned int i = 0; i < HitAndMissMC->numEvents();
         ++i) { // loop over data
      Event event(HitAndMissMC->event(i));

      EventEfficiency = 1.0;
      if (CorrectForEfficiency)
        EventEfficiency = event.efficiency();
      if (EventEfficiency == 0.0) {
        LOG(ERROR) << "DalitzPlot::Fill() | Loop over "
                      "data sample: An event with zero efficiency was found! "
                      "This should not happen! We skip it!";
        continue;
      }
      EventWeight = event.weight();

      /* Construct DataPoint from Event to check if dataPoint is within the
       * phase space boundaries */
      DataPoint point;
      try {
        Kinematics->convert(event, point);
      } catch (BeyondPhsp &ex) { // event outside phase, remove
        continue;
      }

      // Fill branch references with dataPoint
      for (unsigned int j = 0; j < DataPointValues.size(); ++j) {
        DataPointValues[j] = point.value(j);
      }

      tree->Fill();
      bar.next();
    }
    tree->Write();
    delete (tree);
  }
  // write integrals to file
  TVector3 Integral(0, 0, 0);
  Integral.SetX(DataIntegral);
  Integral.Write("data_integral");
  Integral.SetX(PhspIntegral);
  Integral.Write("weighted_phsp_integral");

  tf->Close();
=======
void RootPlotData::write(std::string treePrefix, std::string fileName, std::string option) {

	TFile *tf = new TFile(fileName.c_str(), option.c_str());

	auto varNames = Kinematics->variableNames();
	double EventWeight(1.0);
	double EventEfficiency(1.0);

	double DataIntegral(0.0);
	//===== write data
	if (Data && Data->numEvents() > 0) {
		TTree *tree = new TTree((treePrefix + "_data").c_str(), "DataSample");
		auto DataPointValues = std::vector<double>(varNames.size(), 0.0);
		for (unsigned int i = 0; i < varNames.size(); ++i)
			tree->Branch(varNames.at(i).c_str(), &DataPointValues.at(i), (varNames.at(i) + "/D").c_str());
		tree->Branch("event_weight", &EventWeight, "event_weight/D");
		tree->Branch("event_efficiency", &EventEfficiency, "event_efficiency/D");

		ComPWA::ProgressBar bar(Data->numEvents());
		for (unsigned int i = 0; i < Data->numEvents(); ++i) {  // loop over data
			Event event(Data->event(i));

			EventEfficiency = 1.0;
			if (CorrectForEfficiency) EventEfficiency = event.efficiency();
			if (EventEfficiency == 0.0) {
				LOG(error)<< "DalitzPlot::Fill() | Loop over "
				"data sample: An event with zero efficiency was found! "
				"This should not happen! We skip it!";
				continue;
			}
			EventWeight = event.weight();

			/* Construct DataPoint from Event to check if dataPoint is within the
			 * phase space boundaries */
			DataPoint point;
			try {
				Kinematics->convert(event, point);
			}
			catch (BeyondPhsp &ex) {  // event outside phase, remove
				continue;
			}

			// Fill branch references with dataPoint
			for (unsigned int j = 0; j < DataPointValues.size(); ++j) {
				DataPointValues[j] = point.value(j);
			}

			DataIntegral += point.weight();
			tree->Fill();
			bar.next();
		}
		tree->Write();
		delete (tree);
	}

	double PhspIntegral(0.0);
	//===== write amplitude weights
	if (WeightedPhspMC && WeightedPhspMC->numEvents() > 0) {
		TTree *tree = new TTree((treePrefix + "_weighted_phsp_MC").c_str(), "WeightedPhspMCSample");
		auto DataPointValues = std::vector<double>(varNames.size(), 0.0);
		for (unsigned int i = 0; i < varNames.size(); ++i)
			tree->Branch(varNames.at(i).c_str(), &DataPointValues.at(i), (varNames.at(i) + "/D").c_str());
		tree->Branch("event_weight", &EventWeight, "event_weight/D");
		tree->Branch("event_efficiency", &EventEfficiency, "event_efficiency/D");

		double IntensityWeight(0.0);
		tree->Branch("intensity", &IntensityWeight, "intensity/D");
		std::vector<double> AmplitudeComponentWeights = std::vector<double>(AmplitudeComponents.size(), 0.0);
		unsigned int counter(0);
		for (auto const& amp : AmplitudeComponents) {
			tree->Branch(amp.first.c_str(), &AmplitudeComponentWeights.at(counter),
			    (amp.first + "/D").c_str());
			++counter;
		}

		ComPWA::ProgressBar bar(WeightedPhspMC->numEvents());
		for (unsigned int i = 0; i < WeightedPhspMC->numEvents(); i++) {  // loop over data
			Event event(WeightedPhspMC->event(i));

			EventEfficiency = 1.0;
			if (CorrectForEfficiency) EventEfficiency = event.efficiency();
			if (EventEfficiency == 0.0) {
				LOG(error)<< "DalitzPlot::Fill() | Loop over "
				"data sample: An event with zero efficiency was found! "
				"This should not happen! We skip it!";
				continue;
			}
			EventWeight = event.weight() ;
			PhspIntegral += EventWeight;
			/* Construct DataPoint from Event to check if dataPoint is within the
			 * phase space boundaries */
			DataPoint point;
			try {
				Kinematics->convert(event, point);
			}
			catch (BeyondPhsp &ex) {  // event outside phase, remove
				continue;
			}

			// Fill branch references with dataPoint
			for (unsigned int j = 0; j < DataPointValues.size(); ++j) {
				DataPointValues[j] = point.value(j);
			}

			IntensityWeight = Intensity->intensity(point);
			// Loop over all components that we want to plot
			counter = 0;
			for (auto const& amp : AmplitudeComponents) {
				AmplitudeComponentWeights[counter] = amp.second->intensity(point);
				++counter;
			}

			tree->Fill();
			bar.next();
		}
		tree->Write();
		delete (tree);
	}

	//===== write data
	if (HitAndMissMC && HitAndMissMC->numEvents() > 0) {
		TTree *tree = new TTree((treePrefix + "_hitmiss_MC").c_str(), "Hit&MissMCSample");
		auto DataPointValues = std::vector<double>(varNames.size(), 0.0);
		for (unsigned int i = 0; i < varNames.size(); ++i)
			tree->Branch(varNames.at(i).c_str(), &DataPointValues.at(i), (varNames.at(i) + "/D").c_str());
		tree->Branch("event_weight", &EventWeight, "event_weight/D");
		tree->Branch("event_efficiency", &EventEfficiency, "event_efficiency/D");

		ComPWA::ProgressBar bar(HitAndMissMC->numEvents());
		for (unsigned int i = 0; i < HitAndMissMC->numEvents(); ++i) {  // loop over data
			Event event(HitAndMissMC->event(i));

			EventEfficiency = 1.0;
			if (CorrectForEfficiency) EventEfficiency = event.efficiency();
			if (EventEfficiency == 0.0) {
				LOG(error)<< "DalitzPlot::Fill() | Loop over "
				"data sample: An event with zero efficiency was found! "
				"This should not happen! We skip it!";
				continue;
			}
			EventWeight = event.weight();

			/* Construct DataPoint from Event to check if dataPoint is within the
			 * phase space boundaries */
			DataPoint point;
			try {
				Kinematics->convert(event, point);
			}
			catch (BeyondPhsp &ex) {  // event outside phase, remove
				continue;
			}

			// Fill branch references with dataPoint
			for (unsigned int j = 0; j < DataPointValues.size(); ++j) {
				DataPointValues[j] = point.value(j);
			}

			tree->Fill();
			bar.next();
		}
		tree->Write();
		delete (tree);
	}
	// write integrals to file
	TVector3 Integral(0,0,0);
	Integral.SetX(DataIntegral);
	Integral.Write("data_integral");
	Integral.SetX(PhspIntegral);
	Integral.Write("weighted_phsp_integral");

	tf->Close();
>>>>>>> 417d84dd2... - added RootPlotData class based on RootPlot class.
}

}
}
}
