/**
 * @file GStartIndividual.hpp
 */

/*
 * Copyright (C) Gemfony scientific UG (haftungsbeschraenkt)
 *
 * See the AUTHORS file in the top-level directory for a list of authors.
 *
 * Contact: contact [at] gemfony (dot) com
 *
 * This file is part of the Geneva library collection.
 *
 * Geneva was developed with kind support from Karlsruhe Institute of
 * Technology (KIT) and Steinbuch Centre for Computing (SCC). Further
 * information about KIT and SCC can be found at http://www.kit.edu/english
 * and http://scc.kit.edu .
 *
 * Geneva is free software: you can redistribute and/or modify it under
 * the terms of version 3 of the GNU Affero General Public License
 * as published by the Free Software Foundation.
 *
 * Geneva is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with the Geneva library. If not, see <http://www.gnu.org/licenses/>.
 *
 * For further information on Gemfony scientific and Geneva, visit
 * http://www.gemfony.com .
 */


// Standard header files go here
#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <map>

// Includes check for correct Boost version(s)
#include "common/GGlobalDefines.hpp"

#ifndef GPARABOLOIDINDIVIDUAL2D_HPP_
#define GPARABOLOIDINDIVIDUAL2D_HPP_

// For Microsoft-compatible compilers
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

// ComPWA header files go here
#include <Core/ParameterList.hpp>
#include <Core/Estimator.hpp>

// Geneva header files go here
#include <geneva/GParameterSet.hpp>
#include <geneva/GConstrainedDoubleObject.hpp>

namespace Gem
{
namespace Geneva
{

/*struct namingMapComparator {
  bool operator() (const boost::shared_ptr<GConstrainedDoubleObject> lhs, const boost::shared_ptr<GConstrainedDoubleObject> rhs) const
  {
    //std::cout << "Comperator: " << lhs.get() << "\t" << rhs.get() << std::endl;
    return lhs<rhs;
  }
};*/

/******************************************************************/
/**
 * This individual searches for the minimum of a 2-dimensional parabola.
 * It is part of an introductory example, used in the Geneva manual.
 */
class GStartIndividual :public GParameterSet
{
public:
	/** @brief The default constructor */
	GStartIndividual(std::shared_ptr<ComPWA::IEstimator> data,
			ComPWA::ParameterList list);

	/** @brief A standard copy constructor */
	GStartIndividual(const GStartIndividual&);
	/** @brief The standard destructor */
	virtual ~GStartIndividual();

	bool getPar(ComPWA::ParameterList& val);


	/** @brief A standard assignment operator */
	const GStartIndividual& operator=(const GStartIndividual&);

protected:
	ComPWA::ParameterList parList;
	std::vector<std::string > parNames;

	void updatePar();

	/** @brief Loads the data of another GStartIndividual */
	virtual void load_(const GObject*);
	/** @brief Creates a deep clone of this object */
	virtual GObject* clone_() const;
	/** @brief Loads static data */
	virtual void loadConstantData(boost::shared_ptr<GStartIndividual>);

	/** @brief The actual fitness calculation takes place here. */
	virtual double fitnessCalculation();

private:
	/********************************************************************************************/
	/**
	 * The default constructor. Intentionally private and empty, as it is only needed for
	 * serialization purposes.
	 */
	GStartIndividual();

	/********************************************************************************************/
	// You can add other variables here. Do not forget to serialize them if necessary
	// int myVar;
	std::shared_ptr<ComPWA::IEstimator> theData;

	/** @brief Make the class accessible to Boost.Serialization */
	friend class boost::serialization::access;

	/**************************************************************/
	/**
	 * This function triggers serialization of this class and its
	 * base classes.
	 */
	template<typename Archive>
	void serialize(Archive & ar, const unsigned int) {
		using boost::serialization::make_nvp;
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(GParameterSet);
		// Add other variables here like this:
		// ar & BOOST_SERIALIZATION_NVP(sampleVariable);
	}
	/**************************************************************/

};

/******************************************************************/

} /* namespace Geneva */
} /* namespace Gem */

BOOST_CLASS_EXPORT_KEY(Gem::Geneva::GStartIndividual)

#endif /* GPARABOLOIDINDIVIDUAL2D_HPP_ */
