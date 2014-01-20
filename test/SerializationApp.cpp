/*
 * SerializationTest.cpp
 *
 *  Created on: Jan 17, 2014
 *      Author: weidenka
 */
#include "Core/Parameter.hpp"
#include <boost/serialization/export.hpp>
class MyParent {
public:
	MyParent() {};
	virtual ~MyParent() {};
	int getID(){ return ID;  } ;
	int setID(int num){ ID = num;  } ;

private:
	int ID;
	friend class boost::serialization::access;
	template<class Archive> void serialize(Archive& ar,
			const unsigned int version) {
		ar & BOOST_SERIALIZATION_NVP(ID);
	}
};

class MyChild : public MyParent {
public:
	MyChild() {};
	~MyChild() {};

private:
	friend class boost::serialization::access;
	template<class Archive> void serialize(Archive& ar,
			const unsigned int version) {
		ar & boost::serialization::make_nvp("MyParent", boost::serialization::base_object<MyParent>(*this));
	}
};


template <class T> class classOne
{
public:
	classOne(){};
	classOne(T val) : st(val){};
	virtual ~classOne() {};
	virtual double GetError() =0;
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
BOOST_SERIALIZATION_ASSUME_ABSTRACT(classOne<double>);

template <class T> class classTwo : public classOne<T>
{
public:
	classTwo() : classOne<T>(1.0), error(0){};
	classTwo(T val) : classOne<T>(val), error(val){};
	virtual T GetError() {return error;};

private:
	T error;
	friend class boost::serialization::access;
	template<class archive>
	void serialize(archive& ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(classOne<T>);
		ar & BOOST_SERIALIZATION_NVP(error);
	}

};
//BOOST_CLASS_IMPLEMENTATION(classTwo<double>,boost::serialization::object_serializable);
//BOOST_CLASS_TRACKING(classTwo<double>, boost::serialization::track_always)

int main(int argc, char** argv){
	//	BOOST_SERIALIZATION_ASSUME_ABSTRACT(AbsParameter);
	//	DoubleParameter par("test",2.5,0.0,5.0,0.5);
	//	DoubleParameter parIn();
	//	//write strategy settings
	//	std::ofstream ofs("test.xml");
	//	boost::archive::xml_oarchive oa(ofs);
	//	oa << BOOST_SERIALIZATION_NVP(par);
	//	ofs.close();

	//	std::ifstream ifs("test.xml");
	//	boost::archive::xml_iarchive ia(ifs);
	//	ia >> BOOST_SERIALIZATION_NVP(parIn);
	//	std::cout<<parIn<<std::endl;
	//	ifs.close();

//		SymError<double> err(2.1);
//		SymError<double> errIn(1.0);
//		std::ofstream ofs2("test.xml");
//		boost::archive::xml_oarchive oa2(ofs2,boost::archive::no_header);
//		oa2 << BOOST_SERIALIZATION_NVP(err);
//		ofs2.close();
//		std::ifstream ifs2("test.xml");
//		boost::archive::xml_iarchive ia2(ifs2);
//		ia2 >> BOOST_SERIALIZATION_NVP(errIn);
//		ifs2.close();
	classTwo<double> ch2(0.5);
	std::ofstream ofs2("test2.xml");
	boost::archive::xml_oarchive oa2(ofs2);
	oa2 << BOOST_SERIALIZATION_NVP(ch2);
	ofs2.close();
	std::ifstream ifs2("test2.xml");
	boost::archive::xml_iarchive ia2(ifs2);
	ia2 >> BOOST_SERIALIZATION_NVP(ch2);
	ifs2.close();

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
}

