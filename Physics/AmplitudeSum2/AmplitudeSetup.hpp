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

struct Resonance
{
	std::string m_name;
	//adding type of resonance e.g. couppled channel BW
	std::string m_type;
	std::string m_reference;
	double m_mass; //TODO: struct for borders?
	double m_mass_min;
	double m_mass_max;
	double m_breakup_mom;
	double m_width;
	double m_width_min;
	double m_width_max;
	//adding additional parameters for more complex PDF's, e.g. coupled channel BW
	double m_par1;
	double m_par2;

	double m_strength;
	double m_phase;
	unsigned int m_spin;
	unsigned int m_m;
	unsigned int m_n;

	unsigned int m_daugtherA; //TODO: better reference
	unsigned int m_daugtherB; //TODO: better reference
};

class AmplitudeSetup
{
private:
	std::string m_file;          // ini filename
	std::vector<Resonance> m_resonances;          // resonances
public:
	AmplitudeSetup(const std::string &filename){load(filename);};
	AmplitudeSetup(const AmplitudeSetup& other) : m_file(other.m_file), m_resonances(other.m_resonances) {};
	void load(const std::string &filename);
	void save(const std::string &filename);
	~AmplitudeSetup(){};

	inline const std::string & getFileName() const {return m_file;};
	inline std::vector<Resonance> & getResonances() { return m_resonances; };
};

// Loads amplitude_setup structure from the specified XML file
void AmplitudeSetup::load(const std::string &filename)
{
	// Create an empty property tree object
	using boost::property_tree::ptree;
	ptree pt;

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
			f.m_type = v.second.get<std::string>("type");
			f.m_reference= v.second.get<std::string>("reference");
			f.m_name = v.second.get<std::string>("name");
			f.m_mass = v.second.get<double>("mass");
			f.m_mass_min = v.second.get<double>("mass_min");;
			f.m_mass_max = v.second.get<double>("mass_max");;
			f.m_breakup_mom = v.second.get<double>("breakup_mom");;
			f.m_width = v.second.get<double>("width");;
			f.m_width_min = v.second.get<double>("width_min");;
			f.m_width_max = v.second.get<double>("width_max");;
			f.m_par1 = v.second.get<double>("par1");;
			f.m_par2 = v.second.get<double>("par2");;
			f.m_strength = v.second.get<double>("strength");;
			f.m_phase = v.second.get<double>("phase");;
			f.m_spin = v.second.get<unsigned>("spin");
			f.m_m = v.second.get<unsigned>("m");
			f.m_n = v.second.get<unsigned>("n");
			f.m_daugtherA = v.second.get<unsigned>("daugtherA");
			f.m_daugtherB = v.second.get<unsigned>("daugtherB");
			m_resonances.push_back(f);
		}
	}
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
	BOOST_FOREACH( Resonance const& v, m_resonances ) {
		//     ptree & node = pt.add("amplitude_setup.resonance", "");
		pt.add("amplitude_setup.resonance", "");
		pt.put("name", v.m_name);
		pt.put("type", v.m_type);
		pt.put("mass", v.m_mass);
		pt.put("mass_min", v.m_mass_min);
		pt.put("mass_max", v.m_mass_max);
		pt.put("breakup_mom", v.m_breakup_mom);
		pt.put("width", v.m_width);
		pt.put("width_min", v.m_width_min);
		pt.put("width_max", v.m_width_max);
		pt.put("strength", v.m_strength);
		pt.put("phase", v.m_phase);
		pt.put("spin", v.m_spin);
		pt.put("m", v.m_m);
		pt.put("n", v.m_n);
		pt.put("daughterA", v.m_daugtherA);
		pt.put("daughterB", v.m_daugtherB);
		//if( !v.valid ) node.put("<xmlattr>.invalid", true);
	}

	// Write the property tree to the XML file.
	write_xml(filename, pt);
}

#endif
