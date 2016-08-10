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

#include "GStartIndividual.hpp"

BOOST_CLASS_EXPORT_IMPLEMENT(Gem::Geneva::GStartIndividual)

namespace Gem
{
	namespace Geneva
	{

	using ComPWA::ParameterList;
	using ComPWA::Optimizer::ControlParameter;
	using ComPWA::DoubleParameter;

//	GStartIndividual::GStartIndividual() : GParameterSet(),theData(ControlParameter::Instance())
//	        		{       /* nothing */ }
	GStartIndividual::GStartIndividual() : GParameterSet(),theData()
	{       /* nothing */ }


	/********************************************************************************************/
	/**
	 * The default constructor. This function will add two double parameters to this individual,
	 * each of which has a constrained value range [-10:10].
	 */
	GStartIndividual::GStartIndividual(std::shared_ptr<ControlParameter> data, unsigned int parDim, std::string *name, double* val, double* min, double* max, double* err)
	: GParameterSet(), theData(data)
	{

		// Set up a GConstrainedDoubleObjectCollection

		// Add bounded double objects
		for(std::size_t i=0; i<parDim; i++) {

			boost::shared_ptr<GConstrainedDoubleObject> gbd_ptr(new GConstrainedDoubleObject(val[i], min[i], max[i]) );

			// Create a suitable adaptor (sigma=0.1, sigma-adaption=0.5, min sigma=0, max sigma=0,5)
			boost::shared_ptr<GDoubleGaussAdaptor> gdga_ptr(new GDoubleGaussAdaptor(err[i], 0.5, 0., 3*err[i]));
			gdga_ptr->setAdaptionThreshold(1); // Adaption parameters are modified after each adaption
			gdga_ptr->setAdaptionProbability(0.05); // The likelihood for a parameter to be adapted

			// Register the adaptor with GConstrainedDoubleObject objects
			gbd_ptr->addAdaptor(gdga_ptr);

			// Add a GConstrainedDoubleObject object to the collection
			// gbdc_ptr->push_back(gbd_ptr);
			// gpoc_ptr->push_back(gbd_ptr);
			this->push_back(gbd_ptr);
			//parNames.insert( std::pair<boost::shared_ptr<GConstrainedDoubleObject>,std::string>(gbd_ptr,name[i]) );
			parNames.push_back(name[i]);
		}
	}

  /********************************************************************************************/
  /**
   * A standard copy constructor. All real work is done by the parent class.
   *
   * @param cp A copy of another GStartIndividual
   */
  //GStartIndividual::GStartIndividual(const GStartIndividual& cp)
  //  : GParameterSet(cp),parNames(cp.parNames)
  //  , theData(ControlParameter::Instance())
  //{ /* nothing */ }
  GStartIndividual::GStartIndividual(const GStartIndividual& cp)
	: GParameterSet(cp),parNames(cp.parNames)
	, theData(cp.theData)
	{ /* nothing */ }

	/********************************************************************************************/
	/**
	 * The standard destructor. Note that you do not need to care for the parameter objects
	 * added in the constructor. Upon destruction, they will take care of releasing the allocated
	 * memory.
	 */
	GStartIndividual::~GStartIndividual()
	{ /* nothing */	}

	/********************************************************************************************/
	/**
	 * Access the Parameter
	 *
	 * @param val The ParameterList to fill
	 * @return bool if valid
	 */
	bool GStartIndividual::getPar(ParameterList& val){
		GStartIndividual::conversion_iterator<GConstrainedDoubleObject> it(this->end());
		unsigned int i=0;
		for(it=this->begin(); it!=this->end(); ++it) {
			std::string test;
			if(i<parNames.size())
				test = parNames.at(i);
			val.AddParameter(std::shared_ptr<DoubleParameter>(
					new DoubleParameter(test,(*it)->value())));
			i++;
		}
		return true;
	}

	/********************************************************************************************/
	/**
	 * A standard assignment operator
	 *
	 * @param cp A copy of another GStartIndividual object
	 * @return A constant reference to this object
	 */
	const GStartIndividual& GStartIndividual::operator=(const GStartIndividual& cp){
		GStartIndividual::load_(&cp);
		return *this;
	}

	/********************************************************************************************/
	/**
	 * Loads the data of another GStartIndividual, camouflaged as a GObject.
	 *
	 * @param cp A copy of another GStartIndividual, camouflaged as a GObject
	 */
	void GStartIndividual::load_(const GObject* cp)
	{
		const GStartIndividual *p_load = GObject::gobject_conversion<GStartIndividual>(cp);

		// Load our parent's data
		GParameterSet::load_(cp);

		// No local data
		// sampleVariable = p_load->sampleVariable;
		//theData = ControlParameter::Instance();
	}

	/********************************************************************************************/
	/**
	 * Creates a deep clone of this object
	 *
	 * @return A deep clone of this object, camouflaged as a GObject
	 */
	GObject* GStartIndividual::clone_() const {
		return new GStartIndividual(*this);
	}

	/********************************************************************************************/
	/**
	 * The actual fitness calculation takes place here.
	 *
	 * @return The value of this object
	 */
	double GStartIndividual::fitnessCalculation(){
		//double result = 0.;

		// Extract the GConstrainedDoubleObjectCollection object. In a realistic scenario, you might want
		// to add error checks here upon first invocation.

		GStartIndividual::conversion_iterator<GConstrainedDoubleObject> it(this->end());
		ParameterList minPar;
		unsigned int i=0;
		//std::cout << "STart" <<std::endl;
		for(it=this->begin(); it!=this->end(); ++it) {
			std::string test;
			if(i<parNames.size())
				test = parNames.at(i);
			//std::cout << "ParName: " << test << " Value: "<< (*it)->value() << std::endl;
			//boost::shared_ptr<GConstrainedDoubleObject> ptr = (*it);
			//std::cout << "Value: " << (*it)->value() << " Iterator: "<< (*it).get() << std::endl;

			minPar.AddParameter(std::shared_ptr<DoubleParameter>(
					new DoubleParameter(test,(*it)->value())));
			//  minPar.AddParameter(DoubleParameter("test",(*it)->value(),0,0,0));

			i++;
		}

		if(!theData){
			std::cout << "No Data!" << std::endl; //TODO: exception
			return 0;
		}
		return theData->controlParameter(minPar);
	}

	/********************************************************************************************/
	/**
	 * Loads all static data needed in client mode.
	 */
	void GStartIndividual::loadConstantData(boost::shared_ptr<GStartIndividual>){
		std::cout << "Load data" << std::endl;
		//  theData = ControlParameter::Instance();
	}

	/********************************************************************************************/

	} /* namespace Geneva */
} /* namespace Gem */
