/*
 * AmplitudeSetup.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: weidenka
 */

#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"

namespace ComPWA {
namespace Physics {
namespace AmplitudeSum {

using boost::property_tree::ptree;

AmplitudeSetup::AmplitudeSetup() {
	BOOST_LOG_TRIVIAL(debug) << "AmplitudeSetup::AmplitudeSetup() no filename passed. "
			"Creating AmplitudeSetup with a non resonant component only!";
	basicConf tmp;
	tmp.m_name=="nonRes";
	tmp.m_strength=1.0;
	tmp.m_strength_fix=1.0;
	tmp.m_strength_min=0.0;
	tmp.m_strength_max=2.0;
	tmp.m_phase = 0.0;
	tmp.m_phase_fix=1.0;
	tmp.m_phase_min=-100;
	tmp.m_phase_max=100;
	tmp.m_enable=1;
	m_nonRes.push_back(tmp);
};

AmplitudeSetup::AmplitudeSetup(const std::string &filename) {
	if(filename=="") {
		BOOST_LOG_TRIVIAL(debug) << "AmplitudeSetup::AmplitudeSetup() no filename passed. "
				"Creating AmplitudeSetup with a non resonant component only!";
		basicConf tmp;
		tmp.m_name=="nonRes";
		tmp.m_strength=1.0;
		tmp.m_strength_fix=1.0;
		tmp.m_strength_min=0.0;
		tmp.m_strength_max=2.0;
		tmp.m_phase = 0.0;
		tmp.m_phase_fix=1.0;
		tmp.m_phase_min=-100;
		tmp.m_phase_max=100;
		tmp.m_enable=1;
		m_nonRes.push_back(tmp);
		return;
	}
	load(filename);
};

// Loads amplitude_setup structure from the specified XML file
void AmplitudeSetup::load(const std::string &filename)
{
	// Create an empty property tree object
	using boost::property_tree::ptree;
	//	ptree pt;

	m_filePath=filename;
	// Load the XML file into the property tree. If reading fails
	// (cannot open file, parse error), an exception is thrown.
	read_xml(filename, pt, boost::property_tree::xml_parser::trim_whitespace);
	// Get the filename and store it in the m_file variable.
	// Note that we construct the path to the value by separating
	// the individual keys with dots. If dots appear in the keys,
	// a path type with a different separator can be used.
	// If the amplitude_setup.filename key is not found, an
	// exception is thrown.
	m_file = pt.get<std::string>("amplitude_setup.filename");

	// Iterate over the amplitude_setup.resonances section and
	// store all found resonance names in the m_name set.
	// The get_child() function returns a reference to the child
	// at the specified path; if there is no such child, it throws.
	// Property tree iterators are models of BidirectionalIterator.
	BOOST_FOREACH( ptree::value_type const& v, pt.get_child("amplitude_setup") ) {
		if( v.first == "BreitWigner" ) {
			m_breitWigner.push_back(BreitWignerConf(v.second));
		} else if( v.first == "Flatte" ) {
			m_flatte.push_back(FlatteConf(v.second));
		} else if( v.first == "Flatte3Ch" ) {
			m_flatte3ch.push_back(Flatte3ChConf(v.second));
		} else if( v.first == "nonRes" ) {
			m_nonRes.push_back(basicConf(v.second));
		} else if( v.first == "filename") {
		} else
			throw std::runtime_error("AmpltiudeSetup::load() unknown type of resonance: "+v.first);
	}
	BOOST_LOG_TRIVIAL(info) << "AmplitudeSetup::load() file " << filename
			<< " with " << getNres() <<" resonances!";
	return;
}

// Saves the debug_settings structure to the specified XML file
void AmplitudeSetup::save(const std::string &filename)
{
	// Put log filename in property tree
	pt.put("amplitude_setup.filename", filename);

	// Iterate over the modules in the set and put them in the
	// property tree. Note that the put function places the new
	// key at the end of the list of keys. This is fine most of
	// the time. If you want to place an item at some other place
	// (i.e. at the front or somewhere in the middle), this can
	// be achieved using a combination of the insert and put_own
	// functions.
	BOOST_LOG_TRIVIAL(debug) << "AmplitudeSetup: Saving resonances to "<<filename;
	BOOST_FOREACH( ptree::value_type &v, pt.get_child("amplitude_setup") ) {
		if( v.first == "BreitWigner" ) {
			std::string name = v.second.get<std::string>("name");
			int id = -1;
			for(unsigned int i=0; i<m_breitWigner.size(); i++){
				if(m_breitWigner.at(i).m_name==name) id=i;
			}
			if(id==-1) continue;
			m_breitWigner.at(id).put(v.second);
		}
		if( v.first == "Flatte" ) {
			std::string name = v.second.get<std::string>("name");
			int id = -1;
			for(unsigned int i=0; i<m_flatte.size(); i++){
				if(m_flatte.at(i).m_name==name) id=i;
			}
			if(id==-1) continue;
			m_flatte.at(id).put(v.second);
		}
		if( v.first == "Flatte3Ch" ) {
			std::string name = v.second.get<std::string>("name");
			int id = -1;
			for(unsigned int i=0; i<m_flatte3ch.size(); i++){
				if(m_flatte3ch.at(i).m_name==name) id=i;
			}
			if(id==-1) continue;
			m_flatte3ch.at(id).put(v.second);
		}
		if( v.first == "nonRes" ) {
			std::string name = v.second.get<std::string>("name");
			int id = -1;
			for(unsigned int i=0; i<m_nonRes.size(); i++){
				if(m_nonRes.at(i).m_name==name) id=i;
			}
			if(id==-1) continue;
			m_nonRes.at(id).put(v.second);
		}
	}
	boost::property_tree::xml_writer_settings<char> settings('\t', 1);//new line at the end
	// Write the property tree to the XML file.
	write_xml(filename, pt,std::locale(), settings);
	return;
}

} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */
