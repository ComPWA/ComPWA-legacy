//-------------------------------------------------------------------------------
// Copyright (c) 2013 michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Reader for data in Root-files.
/*! \class RootReader
 * @file RootReader.hpp
 * This class reads event-based data from root-files. It implements the
 * interface of Data.hpp.
 */

#ifndef _RootReader_HPP
#define _RootReader_HPP

#include <vector>
#include <map>
#include <memory>
#include <string>
#include <utility>

//PWA-Headers
#include "DataReader/Data.hpp"
#include "Core/Event.hpp"
#include "Core/DataPoint.hpp"

//Root-Headers
#include "TMath.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TROOT.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TRandom3.h"

class RootReader : public Data {

public:
	//! Default Constructor: Empty object, can be filled
	RootReader();
	RootReader(TTree* tr, const bool binned);
	/*! Read constructor.
	 *
	 * @param inRootFile Input file name
	 * @param inTreeName Name of tree in input file
	 * @param binned	 Create binning
	 */
	RootReader(const std::string inRootFile, const std::string inTreeName="data", const bool binned=false);

	virtual const std::vector<std::string>& getVariableNames();

//  virtual void writeToFile();
    virtual void pushEvent(const Event& evt) {fEvents.push_back(evt);};
    virtual Event& getEvent(const int);
    virtual allMasses getMasses(const unsigned int startEvent=0, unsigned int nEvents=0);
    virtual const int getBin(const int, double&, double&);
    //virtual const int getEvent(const int, TLorentzVector& , TLorentzVector& , double&);
    virtual void writeData(std::string file="",std::string trName="");
    virtual void Clear();

	virtual const unsigned int getNEvents() const {return fEvents.size();};
	virtual const unsigned int getNBins() const {return fmaxBins;};

	virtual std::vector<Event> getEvents() {return fEvents; }
	virtual void Add(Data& otherSample){
		std::vector<Event> otherEvents = otherSample.getEvents();
		fEvents.insert(fEvents.end(), otherEvents.begin(), otherEvents.end());
	}
	/** Destructor */
	virtual ~RootReader();

	//! select only first @param newSize events from full sample
	virtual void reduce(unsigned int newSize){
		if(newSize >= fEvents.size()) {
			BOOST_LOG_TRIVIAL(error) << "RooReader::reduce() requested size too large, cant reduce sample!";
			return;
		}
		fEvents.resize(newSize);
	}
	std::shared_ptr<Data> rndSubSet(unsigned int size, std::shared_ptr<Generator> gen);

	//! Set efficiency value for all stored events. Efficiency is taken from Efficiency object.
	void setEfficiency(std::shared_ptr<Efficiency> eff);
	//! Reset effciencies of all events
	void resetEfficiency(double e=1.);
	//! Reset weights
	void resetWeights(double w=1.);
	//! Weights set?
	bool hasWeights();

protected:
	void read();
	bool _readFlag;

	//input tree
	std::string fileName;
	std::string treeName;
	TFile* fFile;
	TTree* fTree;
	TClonesArray* fParticles;
	double feventWeight;
	double feventEff;
	int fCharge;
	int fFlavour;

	//binning
	bool fBinned;
	unsigned int fmaxBins;
	std::map<int, std::pair<double,double> > fBins;

	//storage of events
	std::vector<std::string> fVarNames;
	std::vector<Event> fEvents;
//	unsigned int fmaxEvents;
	//  unsigned int fEvent;

	virtual void storeEvents();
	virtual void bin();

};

#endif /* _RootReader_HPP */
