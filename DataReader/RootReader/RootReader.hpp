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
#include "Core/Event.hpp"
#include "Core/DataPoint.hpp"
#include "DataReader/Data.hpp"
#include "DataReader/DataCorrection.hpp"

//Root-Headers
//#include "TMath.h"
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
	RootReader(TTree* tr, int size=-1 , const bool binned=false);
	/*! Read constructor.
	 *
	 * @param inRootFile Input file name
	 * @param inTreeName Name of tree in input file
	 * @param size 		Number of events to read in
	 * @param binned	 Create binning
	 */
	RootReader(const std::string inRootFile, const std::string inTreeName="data",
			int size=-1, const bool binned=false);

	virtual RootReader* Clone() const;
	virtual RootReader* EmptyClone() const;
	virtual const std::vector<std::string>& getVariableNames();

    virtual void pushEvent(const Event& evt);
    virtual Event& getEvent(const int);
    virtual allMasses getMasses(const unsigned int startEvent=0, unsigned int nEvents=0);
	virtual ParameterList& getListOfData();

    virtual const int getBin(const int, double&, double&);
    virtual void writeData(std::string file="",std::string trName="");
    virtual void Clear();

    virtual void applyCorrection(DataCorrection& corr);
    //! Number of events
	virtual const unsigned int getNEvents() const {return fEvents.size();};
    //! Number of bins
	virtual const unsigned int getNBins() const {return fmaxBins;};
    //! Get sample as vector of events
	virtual std::vector<Event> getEvents() {return fEvents; }
	//! Get sample as vector of dataPoints
	virtual std::vector<dataPoint> getDataPoints();
	//! Append other sample
	virtual void Add(Data& otherSample);
	//! Destructor
	virtual ~RootReader();
	//! Remove all events outside PHSP
	virtual void reduceToPhsp();
	//! Select only first @param newSize events from full sample
	virtual void reduce(unsigned int newSize);
	//! Select random subset of events
//	std::shared_ptr<Data> rndSubSet(unsigned int size, std::shared_ptr<Generator> gen);

	//! Add resolution to all stored events
	void setResolution(std::shared_ptr<Resolution> res);

	//! Set efficiency value for all stored events. Efficiency is taken from Efficiency object.
	void setEfficiency(std::shared_ptr<Efficiency> eff);
	//! Reset effciencies of all events
	void resetEfficiency(double e=1.);
	//! Reset weights
	void resetWeights(double w=1.);
	//! Get maximum weight
	double getMaxWeight() const;
	//! Weights set?
	bool hasWeights();

protected:
	void read();
	int readSize;
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
	double maxWeight;

	virtual void storeEvents();
	virtual void bin();

};

#endif /* _RootReader_HPP */
