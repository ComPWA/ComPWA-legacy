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

//	GStartIndividual::GStartIndividual() : GParameterSet(),theData(ControlParameter::Instance())
//	        		{       /* nothing */ }

	GStartIndividual::GStartIndividual() : GParameterSet(),theData()
	{       /* nothing */ }

	/********************************************************************************************/
	/**
	 * The default constructor. This function will add two double parameters to this individual,
	 * each of which has a constrained value range [-10:10].
	 */
	GStartIndividual::GStartIndividual(std::shared_ptr<ComPWA::ControlParameter> data, ComPWA::ParameterList list)
	: GParameterSet(), theData(data), parList(list)
	{
		for(std::size_t i=0; i<parList.GetNDouble(); i++) {
			std::shared_ptr<ComPWA::DoubleParameter> p = parList.GetDoubleParameter(i);
			if(p->IsFixed()) continue;
			double val = p->GetValue();
			double min = -1.79768e+307;
			double max =  1.79768e+307;
			double err = val;
			if(p->HasError()) err = p->GetError();
			if(p->UseBounds()){
				min=p->GetMinValue();
				max=p->GetMaxValue();
			}
			boost::shared_ptr<GConstrainedDoubleObject> gbd_ptr(
					new GConstrainedDoubleObject(val, min, max) );

			// Create a suitable adaptor (sigma=0.1, sigma-adaption=0.5, min sigma=0, max sigma=0,5)
			boost::shared_ptr<GDoubleGaussAdaptor> gdga_ptr(new GDoubleGaussAdaptor(
					err, 0.5, 0., 3*err));
			gdga_ptr->setAdaptionThreshold(1); // Adaption parameters are modified after each adaption
			gdga_ptr->setAdaptionProbability(0.05); // The likelihood for a parameter to be adapted

			// Register the adaptor with GConstrainedDoubleObject objects
			gbd_ptr->addAdaptor(gdga_ptr);

			// Add a GConstrainedDoubleObject object to the collection
			// gbdc_ptr->push_back(gbd_ptr);
			// gpoc_ptr->push_back(gbd_ptr);
			this->push_back(gbd_ptr);
			//parNames.insert( std::pair<boost::shared_ptr<GConstrainedDoubleObject>,std::string>(gbd_ptr,name[i]) );
			parNames.push_back(p->GetName());
		}
		BOOST_LOG_TRIVIAL(info) << "GStartIndividual::GStartIndividual() | "
				<<parNames.size()<<" Parameters were added for minimization!";
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
	, theData(cp.theData), parList(cp.parList)
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
	bool GStartIndividual::getPar(ComPWA::ParameterList& val){
		updatePar();
		val = ComPWA::ParameterList(parList);
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
		updatePar();
		return theData->controlParameter(parList);
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

	void GStartIndividual::updatePar(){
		BOOST_LOG_TRIVIAL(debug) << "GStartIndividual::updatePar() | setting new parameters!";
		GStartIndividual::conversion_iterator<GConstrainedDoubleObject> it(this->end());
		it = this->begin();
		for(unsigned int i=0; i<parList.GetNDouble(); ++i) {
			std::shared_ptr<ComPWA::DoubleParameter> p = parList.GetDoubleParameter(i);
			if(p->IsFixed()) continue;
			p->SetValue( (*it)->value() );
			++it;
		}

	}
	} /* namespace Geneva */
} /* namespace Gem */
