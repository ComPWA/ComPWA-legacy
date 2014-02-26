/*
 * AmplitudeSetup.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: weidenka
 */




#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"

// Loads amplitude_setup structure from the specified XML file
void AmplitudeSetup::load(const std::string &filename)
{
	// Create an empty property tree object
	using boost::property_tree::ptree;
	ptree pt;

	m_filePath=filename;
	// Load the XML file into the property tree. If reading fails
	// (cannot open file, parse error), an exception is thrown.
	read_xml(filename, pt);
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
		if( v.first == "resonance" ) {
			Resonance f;
//			f.m_reference= v.second.get<std::string>("reference");
			f.m_enable = v.second.get<bool>("enable");
			f.m_name = v.second.get<std::string>("name");
			f.m_mass = v.second.get<double>("mass");
			f.m_mass_min = v.second.get<double>("mass_min");;
			f.m_mass_max = v.second.get<double>("mass_max");;
			f.m_mesonRadius= v.second.get<double>("mesonRadius");;
			f.m_width = v.second.get<double>("width");;
			f.m_width_min = v.second.get<double>("width_min");;
			f.m_width_max = v.second.get<double>("width_max");;
			f.m_norm= v.second.get<double>("norm");;
			f.m_strength = v.second.get<double>("strength");;
			f.m_phase = v.second.get<double>("phase");;
			f.m_spin = v.second.get<unsigned>("spin");
			f.m_m = v.second.get<unsigned>("m");
			f.m_n = v.second.get<unsigned>("n");
			f.m_daugtherA = v.second.get<unsigned>("daugtherA");
			f.m_daugtherB = v.second.get<unsigned>("daugtherB");
			m_resonances.push_back(f);
			if(f.m_enable) nResEnabled++;
			nRes++;
		}
		if( v.first == "resonanceFlatte" ) {
			ResonanceFlatte f;
//			f.m_reference= v.second.get<std::string>("reference");
			f.m_enable = v.second.get<bool>("enable");
			f.m_name = v.second.get<std::string>("name");
			f.m_mass = v.second.get<double>("mass");
			f.m_mass_min = v.second.get<double>("mass_min");;
			f.m_mass_max = v.second.get<double>("mass_max");;
			f.m_mesonRadius= v.second.get<double>("mesonRadius");;
			f.m_width = v.second.get<double>("width");;
			f.m_width_min = v.second.get<double>("width_min");;
			f.m_width_max = v.second.get<double>("width_max");;
			f.m_norm= v.second.get<double>("norm");;
			f.m_strength = v.second.get<double>("strength");;
			f.m_phase = v.second.get<double>("phase");;
			f.m_spin = v.second.get<unsigned>("spin");
			f.m_m = v.second.get<unsigned>("m");
			f.m_n = v.second.get<unsigned>("n");
			f.m_daugtherA = v.second.get<unsigned>("daugtherA");
			f.m_daugtherB = v.second.get<unsigned>("daugtherB");
			f.m_coupling = v.second.get<double>("coupling");
			f.m_couplingHidden = v.second.get<double>("couplingHidden");
			f.m_hiddenParticle1 = v.second.get<std::string>("hiddenParticle1");
			f.m_hiddenParticle2 = v.second.get<std::string>("hiddenParticle2");
			m_resonancesFlatte.push_back(f);
			if(f.m_enable) nResEnabled++;
			nRes++;
		}
	}
	BOOST_LOG_TRIVIAL(info) << "AmplitudeSetup::load() file " << filename
			<< " with " << nRes	<< "("<<nResEnabled<<") resonances all(enabled)!";
	return;
}

// Saves the debug_settings structure to the specified XML file
void AmplitudeSetup::save(const std::string &filename)
{
	// Create an empty property tree object
	using boost::property_tree::ptree;
	ptree pt;

	// Put log filename in property tree
	pt.put("amplitude_setup.filename", m_file);

	// Iterate over the modules in the set and put them in the
	// property tree. Note that the put function places the new
	// key at the end of the list of keys. This is fine most of
	// the time. If you want to place an item at some other place
	// (i.e. at the front or somewhere in the middle), this can
	// be achieved using a combination of the insert and put_own
	// functions.

	BOOST_LOG_TRIVIAL(debug) << "AmplitudeSetup: Saving resonences "<<m_resonances.size()<<"x BW and "<<m_resonancesFlatte.size()<<"x flatte!";
	BOOST_FOREACH( Resonance const& v, m_resonances ) {
		ptree & node = pt.add("amplitude_setup.resonance", "");
		node.put("enable", v.m_enable);
		node.put("name", v.m_name);
		node.put("mass", v.m_mass);
		node.put("mass_min", v.m_mass_min);
		node.put("mass_max", v.m_mass_max);
		node.put("mesonRadius", v.m_mesonRadius);
		node.put("width", v.m_width);
		node.put("width_min", v.m_width_min);
		node.put("width_max", v.m_width_max);
		node.put("norm",v.m_norm);
		node.put("strength", v.m_strength);
		node.put("phase", v.m_phase);
		node.put("spin", v.m_spin);
		node.put("m", v.m_m);
		node.put("n", v.m_n);
		node.put("daugtherA", v.m_daugtherA);
		node.put("daugtherB", v.m_daugtherB);
//		//if( !v.valid ) node.put("<xmlattr>.invalid", true);
	}
	BOOST_FOREACH( ResonanceFlatte const& v, m_resonancesFlatte ) {
		ptree & node = pt.add("amplitude_setup.resonanceFlatte", "");
		node.put("enable", v.m_enable);
		node.put("name", v.m_name);
		node.put("mass", v.m_mass);
		node.put("mass_min", v.m_mass_min);
		node.put("mass_max", v.m_mass_max);
		node.put("mesonRadius", v.m_mesonRadius);
		node.put("width", v.m_width);
		node.put("width_min", v.m_width_min);
		node.put("width_max", v.m_width_max);
		node.put("norm",v.m_norm);
		node.put("coupling", v.m_coupling);
		node.put("couplingHidden", v.m_couplingHidden);
		node.put("hiddenParticle1", v.m_hiddenParticle1);
		node.put("hiddenParticle2", v.m_hiddenParticle2);
		node.put("strength", v.m_strength);
		node.put("phase", v.m_phase);
		node.put("spin", v.m_spin);
		node.put("m", v.m_m);
		node.put("n", v.m_n);
		node.put("daugtherA", v.m_daugtherA);
		node.put("daugtherB", v.m_daugtherB);
//		//if( !v.valid ) node.put("<xmlattr>.invalid", true);
	}


	// Write the property tree to the XML file.
	write_xml(filename, pt);
}
