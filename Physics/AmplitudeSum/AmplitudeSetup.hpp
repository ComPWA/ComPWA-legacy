//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding parameters for flatte description
//-------------------------------------------------------------------------------
//! XML config parser for Amplitude Setup
/*!
 * @file AmplitudeSetup.hpp
 *\class AmplitudeSetup
 * This class is used to load an Amplitude configuration provided in an XML
 * file. As a Helper the struct Resonance is also provided, of which this class
 * stores a vector.
 */

#ifndef _AMPLITUDESETUP_HPP
#define _AMPLITUDESETUP_HPP

// Standard header files go here
#include <vector>
#include <string>
#include <memory>

// Boost header files go here
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/log/trivial.hpp>

#include "Core/ParameterList.hpp"
#include "Physics/AmplitudeSum/NonResonant.hpp"
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"
#include "Physics/AmplitudeSum/AmpFlatteRes3Ch.hpp"

using namespace boost::log;
using boost::property_tree::ptree;

class AmplitudeSetup
{
private:
	std::string m_file;          // ini filename
	std::string m_filePath;
	std::vector<BreitWignerConf> m_breitWigner;          // resonances
	std::vector<FlatteConf> m_flatte;          // resonances
	std::vector<Flatte3ChConf> m_flatte3ch;          // resonances
	std::vector<basicConf> m_nonRes;
	boost::property_tree::ptree pt;


public:
	AmplitudeSetup();
	AmplitudeSetup(const std::string &filename);
	void load(const std::string &filename);
	void save(const std::string &filename);
	~AmplitudeSetup(){};
	void update(ParameterList par){
		for(std::vector<BreitWignerConf>::iterator reso=getBreitWigner().begin(); reso!=getBreitWigner().end(); reso++)
			(*reso).update(par);
		for(std::vector<FlatteConf>::iterator reso=getFlatte().begin(); reso!=getFlatte().end(); reso++)
			(*reso).update(par);
		for(std::vector<Flatte3ChConf>::iterator reso=getFlatte3Ch().begin(); reso!=getFlatte3Ch().end(); reso++)
			(*reso).update(par);
		for(std::vector<basicConf>::iterator reso=getNonRes().begin(); reso!=getNonRes().end(); reso++)
			(*reso).update(par);
	}

	inline unsigned int getNres() {return m_breitWigner.size()+m_flatte.size()+m_nonRes.size(); };
	inline const std::string & getFileName() const {return m_file;};
	inline const std::string & getFilePath() const {return m_filePath;};
	inline std::vector<BreitWignerConf> & getBreitWigner() { return m_breitWigner; };
	inline std::vector<FlatteConf> & getFlatte() { return m_flatte; };
	inline std::vector<Flatte3ChConf> & getFlatte3Ch() { return m_flatte3ch; };
	inline std::vector<basicConf> & getNonRes() { return m_nonRes; };

};
#endif
