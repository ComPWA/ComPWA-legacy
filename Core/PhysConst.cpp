/*
 * PhysConst.cpp
 *
 *  Created on: Oct 22, 2013
 *      Author: weidenka
 */


#include "Core/PhysConst.hpp"
// Boost header files go here
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

PhysConst* PhysConst::inst = NULL;

PhysConst::PhysConst(){
	//error code
	id.push_back(-999); name.push_back("error"); mass.push_back(-999); width.push_back(-999); J.push_back(999); P.push_back(false); C.push_back(false);
	nameConst.push_back("error"); valueConst.push_back(-999); errorConst.push_back(-999);

	particleFileName = "Physics/DPKinematics/particles.xml";//TODO: dont hardcode datafile
	constantFileName = "Physics/DPKinematics/physConstants.xml";//TODO: dont hardcode datafile

	flag_readFile=1;
	return;
}
void PhysConst::readFile(){
	std::cout<<"PhysConst: reading file with particle information "<<particleFileName<<std::endl;

	// Create an empty property tree object
	using boost::property_tree::ptree;
	ptree pt;
	ptree pt2;

	// Load the XML file into the property tree. If reading fails
	// (cannot open file, parse error), an exception is thrown.
	read_xml(particleFileName, pt);

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
//		std::cout<<"PhysConst DEBUG adding particle: "<<_name<<" mass="<<_mass<<" width="<<_width<<" J=" <<_J<<" P="<<_P<< " C="<<_C<<std::endl;
	}

	read_xml(constantFileName, pt);
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
//		std::cout<<"PhysConst DEBUG adding particle: "<<_name<<" mass="<<_mass<<" width="<<_width<<" J=" <<_J<<" P="<<_P<< " C="<<_C<<std::endl;
	}

	flag_readFile=0;
	return;
}
int PhysConst::findParticle(int idid){
	if(flag_readFile) readFile();
	for(unsigned int i=0; i<name.size(); i++){
		if(id[i]==idid) return i;
	}
	return 0; //error particle not found
}
int PhysConst::findConstant(std::string nnn){
	if(flag_readFile) readFile();
	for(unsigned int i=0; i<nameConst.size(); i++)
		if(nameConst[i]==nnn) return i;
	return 0; //error particle not found
}
int PhysConst::findParticle(std::string nnn){
	if(flag_readFile) readFile();
	for(unsigned int i=0; i<name.size(); i++)
		if(name[i]==nnn) return i;
	return 0; //error particle not found
}
double PhysConst::getMass(std::string nnn)	 	{ return mass[findParticle(nnn)]; } ;
double PhysConst::getWidth(std::string nnn) 	{ return width[findParticle(nnn)]; } ;
unsigned int PhysConst::getJ(std::string nnn) 	{ return J[findParticle(nnn)]; } ;
bool PhysConst::getP(std::string nnn) 			{ return P[findParticle(nnn)]; } ;
bool PhysConst::getC(std::string nnn) 			{ return C[findParticle(nnn)]; } ;

double PhysConst::getMass(int nnn)	 	{ return mass[findParticle(nnn)]; } ;
double PhysConst::getWidth(int nnn) 	{ return width[findParticle(nnn)]; } ;
unsigned int PhysConst::getJ(int nnn) 	{ return J[findParticle(nnn)]; } ;
bool PhysConst::getP(int nnn) 			{ return P[findParticle(nnn)]; } ;
bool PhysConst::getC(int nnn) 			{ return C[findParticle(nnn)]; } ;

double PhysConst::getConstValue(std::string nnn)	 	{ return valueConst[findConstant(nnn)]; } ;
double PhysConst::getConstError(std::string nnn)	 	{ return errorConst[findConstant(nnn)]; } ;
