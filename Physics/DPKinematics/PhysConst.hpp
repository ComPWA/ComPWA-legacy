/*
 * PhysConst.hpp
 *
 *  Created on: Oct 22, 2013
 *      Author: weidenka
 *
 *      PhysConst reads particle properties from an xml file. And provides
 *      it to the framework via a singleton structure. The structure of
 *      the xml file is defined in particleSchema.xsd. It is also intended
 *      to store other physical constants in this way.
 */

#ifndef PHYSCONST_HPP_
#define PHYSCONST_HPP_

#include <iostream>
#include <vector>
#include <Core/Parameter.hpp>

class PhysConst
{
public:
	static PhysConst* instance(){
		if(!inst) inst = new PhysConst();
		return inst;
	}
	double getMass(int);
	double getWidth(int);
	unsigned int getJ(int);
	bool getP(int);
	bool getC(int);
	double getMass(std::string);
	double getWidth(std::string);
	unsigned int getJ(std::string);
	bool getP(std::string);
	bool getC(std::string);
	std::string getName(int);
	int getId(std::string);

	double getConstValue(std::string);
	double getConstError(std::string);

private:
	PhysConst();
	~PhysConst(){};
	PhysConst(PhysConst const&){flag_readFile=1;};
	static PhysConst* inst;
	int findConstant(std::string);
	int findParticle(std::string);
	int findParticle(int);
	void readFile();

	bool flag_readFile;
	std::string pdgFileName;
	std::vector<std::string> name;//TODO: use a particle class to store information
	std::vector<int> id;
	std::vector<double> mass;
	std::vector<double> width;
	std::vector<unsigned int> J;
	std::vector<int> P;
	std::vector<int> C;
	std::vector<std::string> nameConst;
	std::vector<double> valueConst;
	std::vector<double> errorConst;
};


#endif /* PHYSCONST_HPP_ */
