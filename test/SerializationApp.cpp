/*
 * SerializationTest.cpp
 *
 *  Created on: Jan 17, 2014
 *      Author: weidenka
 */
#include "Core/Parameter.hpp"
#include <sstream>
#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/archive/archive_exception.hpp>


template <class T> class classOne
{
public:
	classOne(){};
	classOne(T val) : st(val){};
	virtual ~classOne() {};
	virtual double GetError() {return 1;}
protected:
	T st;
private:
	friend class boost::serialization::access;
	template<class archive>
	void serialize(archive& ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_NVP(st);
	}

};
BOOST_SERIALIZATION_ASSUME_ABSTRACT(classOne);

template <class T> class classTwo : public classOne<T>
{
public:
	virtual ~classTwo() {};
	classTwo(T val) : classOne<T>(val), error(val){};
	virtual T GetError() {return error;};
	classTwo() : classOne<T>(.5), error(1.0){};

private:
	T error;
	friend class boost::serialization::access;
	template<class archive>
	void serialize(archive& ar, const unsigned int version)
	{
		typedef classOne<T> baseClass;
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(baseClass);
		ar & BOOST_SERIALIZATION_NVP(error);
	}

};
//BOOST_CLASS_EXPORT_GUID(classTwo, "classTwo")
//BOOST_CLASS_IMPLEMENTATION(classTwo<double>,boost::serialization::object_serializable);
//BOOST_CLASS_TRACKING(classTwo<double>, boost::serialization::track_always)

int main(int argc, char** argv){
	//===== serilization of DoubleParameter
	DoubleParameter par("test",2.5,0.0,5.0,0.5);
	DoubleParameter parIn("dsfa");
	AsymError<double> parErr(std::pair<double,double>(1.0,2.0));
//	par.SetError(std::shared_ptr<ParError<double>>(new AsymError<double>(std::pair<double,double>(1.0,2.0))));

	std::ofstream ofs2("test.xml");
	boost::archive::xml_oarchive oa2(ofs2);
	oa2 << BOOST_SERIALIZATION_NVP(par);

	std::ifstream ifs("test.xml");
	boost::archive::xml_iarchive ia(ifs);
	ia >> BOOST_SERIALIZATION_NVP(parIn);
	std::cout<<parIn<<std::endl;
	ifs.close();

	//===== serialization of error
//	std::cout<<"SymError: writing"<<std::endl;
//	SymError<double> err(2.1);
//	oa2 << BOOST_SERIALIZATION_NVP(err);
//	std::cout<<"ASymError: writing"<<std::endl;
//	AsymError<double> errA(std::pair<double,double>(1.0,2.0));
//	oa2 << BOOST_SERIALIZATION_NVP(errA);
//
//	ofs2.close();
//	std::ifstream ifs2("test.xml");
//	boost::archive::xml_iarchive ia2(ifs2);
//
//	SymError<double> errIn(1.0);
//	AsymError<double> errAin(std::pair<double,double>(1.0,2.0));
//
//	std::cout<<"SymError: reading"<<std::endl;
//	ia2 >> BOOST_SERIALIZATION_NVP(errIn);
//	std::cout<<errIn.GetError()<<std::endl;
//
//	std::cout<<"ASymError: reading"<<std::endl;
//	ia2 >> BOOST_SERIALIZATION_NVP(errAin);
//	std::cout<<errAin.GetError()<<std::endl;
//	ifs2.close();

	//	classTwo<double> ch2;
	//	std::ofstream ofs3("test2.xml");
	//	boost::archive::xml_oarchive oa3(ofs3);
	//	std::cout<<"Write to file..."<<std::endl;
	//	oa3 << BOOST_SERIALIZATION_NVP(ch2);
	//	ofs3.close();
	//
	//	classTwo<double> chIn;
	//	std::ifstream ifs3("test2.xml");
	//	boost::archive::xml_iarchive ia3(ifs3);
	//	std::cout<<"Read from file..."<<std::endl;
	//	ia3 >> BOOST_SERIALIZATION_NVP(chIn);
	//	ifs3.close();
	//		std::cout<<chIn.GetError()<<std::endl;


	//	MyChild ch2;
	//	ch2.setID(3);
	//	std::ofstream ofs2("test2.xml");
	//	boost::archive::xml_oarchive oa2(ofs2);
	//	oa2 << BOOST_SERIALIZATION_NVP(ch2);
	//	ofs2.close();
	//	std::ifstream ifs2("test2.xml");
	//	boost::archive::xml_iarchive ia2(ifs2);
	//	ia2 >> BOOST_SERIALIZATION_NVP(ch2);
	//	ifs2.close();
	return 0;
}

