//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//	   Peter Weidenkaff - Weights and background fractions
//-------------------------------------------------------------------------------

//! Negative Log Likelihood-Estimator.
/*! \class MinLogLHbkg
 * @file MinLogLHbkg.hpp
 * This class calculates a simple -log(LH) of a intensity and a dataset.
 * Data and Model are provided in the constructor using the Amplitude and Data
 * interfaces. The class itself fulfills the Estimator interface.
 */

#ifndef _MINLOGLHBKG_HPP
#define _MINLOGLHBKG_HPP

#include <vector>
#include <memory>
#include <string>

//PWA-Header
#include "Estimator/Estimator.hpp"
#include "Core/Amplitude.hpp"
#include "DataReader/Data.hpp"
#include "Core/Event.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"

class MinLogLHbkg : public Estimator {

public:
	//! Implementation of ControlParameter::controlParameter
	virtual double controlParameter(ParameterList& minPar);
	/** Create instance of MinLogLHbkg.
	 * A binned efficiency correction is used. We expect that phspSample_ has efficiency values for
	 * each event.
	 *
	 * @param amp_ amplitude
	 * @param data_ data sample
	 * @param phspSample_ phsp sample for normalization. Efficiency values for each point needs
	 *  to be set beforehand.
	 * @param startEvent use #data_ from that position on
	 * @param nEvents number of events to process
	 * @return std::shared_ptr<Data> of existing instance or newly created instance
	 */
	static std::shared_ptr<ControlParameter> createInstance(std::shared_ptr<Amplitude> amp_,
			std::shared_ptr<Data> data_, std::shared_ptr<Data> phspSample_,
			unsigned int startEvent=0, unsigned int nEvents=0);
	/** Create instance of MinLogLHbkg.
	 * An unbinned efficiency correction is applied using #accSample_.
	 *
	 * @param amp_ amplitude
	 * @param data_ data sample
	 * @param phspSample_ phsp sample for normalization
	 * @param accSample_ sample of efficiency applied phsp events for unbinned efficiency correction
	 * @param startEvent use #data_ from that position on
	 * @param nEvents number of events to process
	 * @param sigFrac signal fraction in data sample
	 * @return std::shared_ptr<Data> of existing instance or newly created instance
	 */
	static std::shared_ptr<ControlParameter> createInstance(std::shared_ptr<Amplitude> amp_,
			std::shared_ptr<Data> data_, std::shared_ptr<Data> phspSample_,std::shared_ptr<Data> accSample_,
			unsigned int startEvent=0, unsigned int nEvents=0);
	/** Create instance of MinLogLHbkg.
	 * An unbinned efficiency correction is applied using #accSample_.
	 *
	 * @param amp_ amplitude
	 * @param bkg_ background description amplitude
	 * @param data_ data sample
	 * @param phspSample_ phsp sample for normalization
	 * @param accSample_ sample of efficiency applied phsp events for unbinned efficiency correction
	 * @param startEvent use #data_ from that position on
	 * @param nEvents number of events to process
	 * @param sigFrac signal fraction in data sample
	 * @return std::shared_ptr<Data> of existing instance or newly created instance
	 */
	static std::shared_ptr<ControlParameter> createInstance(std::shared_ptr<Amplitude> amp_,std::shared_ptr<Amplitude> bkg_,
			std::shared_ptr<Data> data_, std::shared_ptr<Data> phspSample_,std::shared_ptr<Data> accSample_,
			unsigned int startEvent, unsigned int nEvents, double sigFrac);

	/** Set new amplitude to existing instance
	 *
	 * @param amp_ amplitude
	 * @param data_ data sample
	 * @param phspSample_ phsp sample for normalization
	 * @param accSample_ sample of efficiency applied phsp events for unbinned efficiency correction
	 * @param startEvent use #data_ from that position on
	 * @param nEvents number of events to process
	 * @param useFuncTr use FunctionTree yes/no?
	 */
	virtual void setAmplitude(std::shared_ptr<Amplitude> amp_, std::shared_ptr<Data> data_,
			std::shared_ptr<Data> phspSample_, std::shared_ptr<Data> accSample_,
			unsigned int startEvent=0, unsigned int nEvents=0, bool useFuncTr=0, double sigFrac=1.);
	/** Set new amplitude to existing instance
	 *
	 * @param amp_ amplitude
	 * @param bkg_ amplitude
	 * @param data_ data sample
	 * @param phspSample_ phsp sample for normalization
	 * @param accSample_ sample of efficiency applied phsp events for unbinned efficiency correction
	 * @param startEvent use #data_ from that position on
	 * @param nEvents number of events to process
	 * @param useFuncTr use FunctionTree yes/no?
	 */
	virtual void setAmplitude(std::shared_ptr<Amplitude> amp_, std::shared_ptr<Amplitude> bkg_,std::shared_ptr<Data> data_,
			std::shared_ptr<Data> phspSample_, std::shared_ptr<Data> accSample_,
			unsigned int startEvent=0, unsigned int nEvents=0, bool useFuncTr=0, double sigFrac=1.);

	//! Check if tree for LH calculation is available
	virtual bool hasTree() { return useFunctionTree; }
	//! Get FunctionTree for LH calculation. Check first if its available using MinLogLHbkg::hasTree().
	virtual std::shared_ptr<FunctionTree> getTree() {
		if(!useFunctionTree)
			throw std::runtime_error("MinLogLHbkg::getTree() you requested a tree which was not constructed!");
		return physicsTree;
	}
	//! Get Amplitude
	virtual std::shared_ptr<Amplitude> getAmplitude() { return amp; }
	//! Destructor
	virtual ~MinLogLHbkg() { };
	//! Should we try to use the function tree? Function tree needs to be implemented in Amplitude
	void setUseFunctionTree(bool t);

protected:
	//! Default Constructor
	MinLogLHbkg() { };
	//! Constructor
	MinLogLHbkg(std::shared_ptr<Amplitude> amp_, std::shared_ptr<Amplitude> bkg_,std::shared_ptr<Data> data_,
			std::shared_ptr<Data> phspSample_, std::shared_ptr<Data> accSample_,
			unsigned int startEvent, unsigned int nEvents, double sigFrac=1.);
	//! Uses ampTree and creates a tree that calculates the full LH
	virtual void iniLHtree();
	//! Sum up all weights in data set
	void calcSumOfWeights();

private:
	//! FunctionTree for Likelihood calculation
	std::shared_ptr<FunctionTree> physicsTree;

	//! Signal decay amplitude
	std::shared_ptr<Amplitude> amp;
	//! FunctionTree calculation the LH normalizaton
	std::shared_ptr<FunctionTree> signalPhspTree;
	//! Signal amplitude tree from #amp with invariant mass from data sample
	std::shared_ptr<FunctionTree> signalTree_amp;
	//! Signal amplitude tree from #amp with invariant masses from phsp/acc sample
	std::shared_ptr<FunctionTree> signalPhspTree_amp;
	//! Background amplitude
	std::shared_ptr<Amplitude> ampBkg;
	//! FunctionTree calculation of background normalization
	std::shared_ptr<FunctionTree> bkgPhspTree;
	//! Background amplitude tree from #ampBkg with invariant mass from data sample
	std::shared_ptr<FunctionTree> bkgTree_amp;
	//! Background amplitude tree from #ampBkg with invariant masses from phsp/acc sample
	std::shared_ptr<FunctionTree> bkgPhspTree_amp;

	//! Data sample
	std::shared_ptr<Data> data;
	//! Phsp sample for normalization
	std::shared_ptr<Data> phspSample;
	//! Phsp with applied efficency for unbinned efficiency correction
	std::shared_ptr<Data> accSample;
	//! Data sample
	allMasses mData;
	//! Phsp sample for normalization
	allMasses mPhspSample;
	//! Phsp with applied efficency for unbinned efficiency correction
	allMasses mAccSample;

	//! Number of events in data sample
	unsigned int nEvts_;
	//! Number of event in phsp sample
	unsigned int nPhsp_;
	//! Process data sample from position #nStartEvt_ on
	unsigned int nStartEvt_;
	//! Number of events to process in data sample
	unsigned int nUseEvt_;
	//! Sum over all weights in data sample
	double sumOfWeights;
	//! Use FunctionTree for LH calculation (yes/no)?
	bool useFunctionTree;
	//! Fraction of signal in data sample
	double signalFraction;

};

#endif /* _MINLOGLHBKG_HPP */
