/*
 * SerializationTest.cpp
 *
 *  Created on: Jan 17, 2014
 *      Author: weidenka
 */
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include <sstream>
#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/archive/archive_exception.hpp>

#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>

using namespace boost::serialization;

using ComPWA::ParameterList;
using ComPWA::DoubleParameter;

int main(int argc, char** argv){

	std::cout<<"Testing boost::serialization"<<std::endl;
	//Serialize a DoubleParameter
	DoubleParameter p1("test",2.5,0.0,5.0,0.5);
	p1.SetError(0.5,1.3);
	DoubleParameter p1In;
	std::ofstream ofs("DoublePar.xml");
	boost::archive::xml_oarchive oa(ofs);
	oa << make_nvp("DoublePar", p1);
	ofs.close();
	std::ifstream ifs("DoublePar.xml");
	boost::archive::xml_iarchive ia(ifs);
	ia >> make_nvp("DoublePar", p1In);
	ifs.close();
	std::cout<<"Serialization of DoubleParameter: value="<<p1In.GetValue()<<std::endl;


	//Serialize a ParameterList
	std::shared_ptr<DoubleParameter> p3(new DoubleParameter("test2",1.5,0.0,5.0,0.5));
	ParameterList list, list2;
	list.AddParameter(p3);
	std::ofstream ofs3("ParList.xml");
	boost::archive::xml_oarchive oa3(ofs3);
	oa3 << boost::serialization::make_nvp("ParList", list);
	ofs3.close();
	std::ifstream ifs3("ParList.xml");
	boost::archive::xml_iarchive ia3(ifs3);
	ia3 >> boost::serialization::make_nvp("ParList", list2);
	std::cout<<"Serialization of ParameterList: NDouble="<<list2.GetNDouble()<<std::endl;

	//Serialize a std::shared_ptr<DoubleParameter>
	std::shared_ptr<DoubleParameter> p2(new DoubleParameter("test2",1.5,0.0,5.0,0.5));
	std::shared_ptr<DoubleParameter> p2In;
	std::ofstream ofs2("ShrDoublePar.xml");
	boost::archive::xml_oarchive oa2(ofs2);
	oa2 << boost::serialization::make_nvp("ShrDoublePar", p2);
	ofs2.close();
	std::ifstream ifs2("ShrDoublePar.xml");
	boost::archive::xml_iarchive ia2(ifs2);
	ia2 >> boost::serialization::make_nvp("ShrDoublePar", p2In);
	std::cout<<"Serialization of std::shared_ptr<DoubleParameter>: name=" << p2In->GetName()<< " value="<<p2In->GetValue()<<std::endl;

	return 0;
}

