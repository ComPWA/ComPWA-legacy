#ifndef OIFDATA_HPP_
#define OIFDATA_HPP_

#include <vector>

class OIFData
{
	///////////////////////////////////////////////////////////////////////
	//friend class boost::serialization::access;

	//template<typename Archive>
	//void serialize(Archive & ar, const unsigned int) {
	//	using boost::serialization::make_nvp;

	//	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(GParameterSet);

		/* Add your own class-variables here in the following way:
			ar & BOOST_SERIALIZATION_NVP(myVar);
			or
			ar & make_nvp("myVar", myVar); // The latter form can be necessary when dealing with templates
		 */
  //}
	///////////////////////////////////////////////////////////////////////

public:

	OIFData()
	  {
	  }

  virtual ~OIFData()
	{ /* nothing */	}

  virtual double controlParameter(const std::vector<double>& minPar) =0;
 
};

#endif
