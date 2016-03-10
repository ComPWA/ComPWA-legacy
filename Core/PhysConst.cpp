/*
 * PhysConst.cpp
 *
 *  Created on: Oct 22, 2013
 *      Author: weidenka
 */

#include <stdlib.h>
#include <sstream>
#include <string>
#include <exception>

#include "Core/PhysConst.hpp"
// Boost header files go here
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/log/trivial.hpp>
using namespace boost::log;

PhysConst* PhysConst::inst = NULL;

PhysConst::PhysConst(){
	//error code
	id.push_back(-999); name.push_back("error"); mass.push_back(-999); width.push_back(-999); J.push_back(999); P.push_back(false); C.push_back(false);
	nameConst.push_back("error"); valueConst.push_back(-999); errorConst.push_back(-999);

	const char* pPath = getenv("COMPWA_DIR");
	std::string path = "";
	try{
		path = std::string(pPath);
	}catch(std::logic_error){
		BOOST_LOG_TRIVIAL(error)<<"Environment Variable COMPWA_DIR not set?"<<std::endl;
	}
	particleFileName = path+"/Physics/particles.xml";
	particleDefaultFileName = path+"/Physics/particlesDefault.xml";
	constantFileName = path+"/Physics/physConstants.xml";
	constantDefaultFileName = path+"/Physics/physDefaultConstants.xml";

	flag_readFile=1;
	return;
}
void PhysConst::readFile(){

	// Create an empty property tree object
	using boost::property_tree::ptree;
	ptree pt;
	ptree pt2;

	// Load the XML file into the property tree. If reading fails
	// (cannot open file, parse error), an exception is thrown.
	//		read_xml(particleFileName, pt);
	//first check if particles.xml existsheck if file exists
	if (FILE *file = std::fopen(particleFileName.c_str(), "r")) {
		fclose(file);
		read_xml(particleFileName, pt);
		BOOST_LOG_TRIVIAL(info) << "PhysConst: reading particle file "<<particleFileName;
		//Otherwise try to load default file
	}else if (FILE *file = std::fopen(particleDefaultFileName.c_str(), "r")) {
		fclose(file);
		read_xml(particleDefaultFileName, pt);
		BOOST_LOG_TRIVIAL(info) << "PhysConst: reading particles default file "<<particleDefaultFileName;
	} else {
		throw std::runtime_error("Could not open default particle file!");
	}

	// Get the filename and store it in the m_file variable.
	// Note that we construct the path to the value by separating
	// the individual keys with dots. If dots appear in the keys,
	// a path type with a different separator can be used.
	// If the amplitude_setup.filename key is not found, an
	// exception is thrown.
	//    m_file = pt.get<std::string>("amplitude_setup.filename");

	// Iterate over the amplitude_setup.resonances section and
	// store all found resonance names in the m_name set.
	// The get_child() function returns a reference to the child
	// at the specified path; if there is no such child, it throws.
	// Property tree iterators are models of BidirectionalIterator.
	int _id;
	std::string _name;
	double _mass;
	double _width;
	unsigned int _J;
	int _P;
	int _C;
	double _error;
	double _value;


	BOOST_FOREACH( ptree::value_type const& v, pt.get_child("particleList") ) {
		_id = -999; _name="error"; _mass=-999; _width=-999; _J=-999; _P=-999; _C=-999;//setting default values
		if( v.first == "particle" || v.first == "particleFlatte") {
			_id= v.second.get<int>("ID",-999);
			_name = v.second.get<std::string>("name","error");
			_mass = v.second.get_child("mass").get<double>("value",-999);
			if( v.second.count("width") != 0) _width = v.second.get_child("width").get<double>("value",-999);//check if node "width" exists
			_J = v.second.get<unsigned int>("J",999);
			_P = v.second.get<int>("P",0);
			_C = v.second.get<int>("C",0);
		}
		if( v.first == "particleFlatte" ) {
			//read parameters which are specific to flatte description here.

		}
		if(_name=="error" && _mass==-999 ) continue;
		id.push_back(_id);
		name.push_back(_name);
		mass.push_back(_mass);
		width.push_back(_width);
		J.push_back(_J);
		P.push_back(_P);
		C.push_back(_C);
		BOOST_LOG_TRIVIAL(debug)<<"PhysConst adding particle: "<<_name<<" mass="<<_mass<<" width="<<_width<<" J=" <<_J<<" P="<<_P<< " C="<<_C;
	}

	//Reading XML file with physics constants
	if (FILE *file = std::fopen(constantFileName.c_str(), "r")) {
		fclose(file);
		read_xml(constantFileName, pt);
		BOOST_LOG_TRIVIAL(info) << "PhysConst: reading file with physical constants"<<constantFileName;
		//Otherwise try to load default file
	}else if (FILE *file = std::fopen(constantDefaultFileName.c_str(), "r")) {
		fclose(file);
		read_xml(constantDefaultFileName, pt);
		BOOST_LOG_TRIVIAL(info) << "PhysConst: reading default file with physical constants"<<constantDefaultFileName;
	} else {
		throw std::runtime_error("Could not open default constants file!");
	}
//	read_xml(constantFileName, pt);
	BOOST_FOREACH( ptree::value_type const& v, pt.get_child("physConstList") ) {
		_name="error"; _error=-999; _value=-999;
		if( v.first == "constant" ){
			_name = v.second.get<std::string>("name","error");
			_value= v.second.get_child("value").get<double>("value",-999);
			_error= v.second.get_child("value").get<double>("error",-999);
		}
		if(_name=="error" && _value==-999 ) continue;
		nameConst.push_back(_name);
		valueConst.push_back(_value);
		errorConst.push_back(_error);
		BOOST_LOG_TRIVIAL(debug)<<"PhysConst adding particle: "<<_name<<" mass="<<_mass<<" width="<<_width<<" J=" <<_J<<" P="<<_P<< " C="<<_C;
	}

	flag_readFile=0;
	return;
}

int PhysConst::findParticle(int idid)
{
	if(flag_readFile) readFile();
	for(unsigned int i=0; i<name.size(); i++){
		if(id[i]==idid) return i;
	}
	return 0; //error particle not found
}

int PhysConst::findConstant(std::string nnn)
{
	if(flag_readFile) readFile();
	for(unsigned int i=0; i<nameConst.size(); i++)
		if(nameConst[i]==nnn) return i;
	return 0; //error particle not found
}
int PhysConst::findParticle(std::string nnn)
{
	if(flag_readFile) readFile();
	for(unsigned int i=0; i<name.size(); i++)
		if(name[i]==nnn) return i;
	throw BadParameter("PhysConst::findParticle() | Particle "+nnn+" not found!");
	return 0;
}

double PhysConst::getMass(std::string nnn)
{
	return mass.at(findParticle(nnn));
}

double PhysConst::getWidth(std::string nnn)
{
	return width.at(findParticle(nnn));
}

unsigned int PhysConst::getJ(std::string nnn)
{
	return J.at(findParticle(nnn));
}

bool PhysConst::getP(std::string nnn)
{
	return P.at(findParticle(nnn));
}

bool PhysConst::getC(std::string nnn)
{
	return C.at(findParticle(nnn));
}

double PhysConst::getMass(int nnn)
{
	return mass.at(findParticle(nnn));
}

double PhysConst::getWidth(int nnn)
{
	return width.at(findParticle(nnn));
}

unsigned int PhysConst::getJ(int nnn)
{
	return J.at(findParticle(nnn));
}

bool PhysConst::getP(int nnn)
{
	return P.at(findParticle(nnn));
}

bool PhysConst::getC(int nnn)
{
	return C.at(findParticle(nnn));
}

double PhysConst::getConstValue(std::string nnn)
{
	return valueConst.at(findConstant(nnn));
}

double PhysConst::getConstError(std::string nnn)
{
	return errorConst.at(findConstant(nnn));
}
