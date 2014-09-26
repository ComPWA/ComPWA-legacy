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
			f.m_mass_fix = v.second.get<bool>("mass_fix");;
			f.m_mass_min = v.second.get<double>("mass_min");;
			f.m_mass_max = v.second.get<double>("mass_max");;
			f.m_width = v.second.get<double>("width");;
			f.m_width_fix = v.second.get<bool>("width_fix");;
			f.m_width_min = v.second.get<double>("width_min");;
			f.m_width_max = v.second.get<double>("width_max");;
			f.m_strength = v.second.get<double>("strength");;
			f.m_strength_fix = v.second.get<bool>("strength_fix");;
			f.m_strength_min = v.second.get<double>("strength_min");;
			f.m_strength_max = v.second.get<double>("strength_max");;
			f.m_phase = v.second.get<double>("phase");;
			f.m_phase_fix = v.second.get<bool>("phase_fix");;
			f.m_phase_min = v.second.get<double>("phase_min");;
			f.m_phase_max = v.second.get<double>("phase_max");;
			f.m_mesonRadius= v.second.get<double>("mesonRadius");;
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
			f.m_mass_fix = v.second.get<bool>("mass_fix");;
			f.m_mass_min = v.second.get<double>("mass_min");;
			f.m_mass_max = v.second.get<double>("mass_max");;
			f.m_strength = v.second.get<double>("strength");;
			f.m_strength_fix = v.second.get<bool>("strength_fix");;
			f.m_strength_min = v.second.get<double>("strength_min");;
			f.m_strength_max = v.second.get<double>("strength_max");;
			f.m_phase = v.second.get<double>("phase");;
			f.m_phase_fix = v.second.get<bool>("phase_fix");;
			f.m_phase_min = v.second.get<double>("phase_min");;
			f.m_phase_max = v.second.get<double>("phase_max");;
			f.m_mesonRadius= v.second.get<double>("mesonRadius");;
			f.m_spin = v.second.get<unsigned>("spin");
			f.m_m = v.second.get<unsigned>("m");
			f.m_n = v.second.get<unsigned>("n");
			f.m_daugtherA = v.second.get<unsigned>("daugtherA");
			f.m_daugtherB = v.second.get<unsigned>("daugtherB");
			f.m_g1 = v.second.get<double>("g1");
			f.m_g1_fix = v.second.get<bool>("g1_fix");;
			f.m_g1_min = v.second.get<double>("g1_min");;
			f.m_g1_max = v.second.get<double>("g1_max");;
			f.m_g2 = v.second.get<double>("g2");
			f.m_g2_part1 = v.second.get<std::string>("g2_part1");
			f.m_g2_part2 = v.second.get<std::string>("g2_part2");
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
		node.put("mass_fix", v.m_mass_fix);
		node.put("mass_min", v.m_mass_min);
		node.put("mass_max", v.m_mass_max);
		node.put("width", v.m_width);
		node.put("width_fix", v.m_width_fix);
		node.put("width_min", v.m_width_min);
		node.put("width_max", v.m_width_max);
		node.put("strength", v.m_strength);
		node.put("strength_fix", v.m_strength_fix);
		node.put("strength_min", v.m_strength_min);
		node.put("strength_max", v.m_strength_max);
		node.put("phase", v.m_phase);
		node.put("phase_fix", v.m_phase_fix);
		node.put("phase_min", v.m_phase_min);
		node.put("phase_max", v.m_phase_max);
		node.put("mesonRadius", v.m_mesonRadius);
		node.put("spin", v.m_spin);
		node.put("m", v.m_m);
		node.put("n", v.m_n);
		node.put("daugtherA", v.m_daugtherA);
		node.put("daugtherB", v.m_daugtherB);
	}
	BOOST_FOREACH( ResonanceFlatte const& v, m_resonancesFlatte ) {
		ptree & node = pt.add("amplitude_setup.resonanceFlatte", "");
		node.put("enable", v.m_enable);
		node.put("name", v.m_name);
		node.put("mass", v.m_mass);
		node.put("mass_fix", v.m_mass_fix);
		node.put("mass_min", v.m_mass_min);
		node.put("mass_max", v.m_mass_max);
		node.put("strength", v.m_strength);
		node.put("strength_fix", v.m_strength_fix);
		node.put("strength_min", v.m_strength_min);
		node.put("strength_max", v.m_strength_max);
		node.put("phase", v.m_phase);
		node.put("phase_fix", v.m_phase_fix);
		node.put("phase_min", v.m_phase_min);
		node.put("phase_max", v.m_phase_max);
		node.put("mesonRadius", v.m_mesonRadius);
		node.put("spin", v.m_spin);
		node.put("m", v.m_m);
		node.put("n", v.m_n);
		node.put("daugtherA", v.m_daugtherA);
		node.put("daugtherB", v.m_daugtherB);
		node.put("g1", v.m_g1);
		node.put("g1_fix", v.m_g1_fix);
		node.put("g1_min", v.m_g1_min);
		node.put("g1_max", v.m_g1_max);
		node.put("g2", v.m_g2);
		node.put("g2_part1", v.m_g2_part1);
		node.put("g2_part2", v.m_g2_part2);

	}


	// Write the property tree to the XML file.
	write_xml(filename, pt);
}
