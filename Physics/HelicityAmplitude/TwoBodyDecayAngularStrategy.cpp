/*
 * TwoBodyDecayAngularStrategy.cpp
 *
 *  Created on: Aug 9, 2016
 *      Author: steve
 */

#include "Core/DataPointStorage.hpp"
#include "Physics/HelicityAmplitude/TwoBodyDecayAngularStrategy.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

TwoBodyDecayAngularStrategy::TwoBodyDecayAngularStrategy(
		std::shared_ptr<TwoBodyDecayAmplitude> tbd_amp, unsigned storage_index) :
    	Strategy(ParType::MCOMPLEX), tbd_amp_(tbd_amp), storage_index_(
    						storage_index)
{

}

TwoBodyDecayAngularStrategy::~TwoBodyDecayAngularStrategy()
{
	// TODO Auto-generated destructor stub
}

//! Pure Virtual interface for streaming info about the strategy
const std::string TwoBodyDecayAngularStrategy::to_str() const
{
	return "WignerD";
}

//! Pure Virtual interface for executing a strategy
bool TwoBodyDecayAngularStrategy::execute(ParameterList& paras,
		std::shared_ptr<AbsParameter>& out)
{

	if ( out->type() != ParType::MCOMPLEX)
		throw BadParameter("TwoBodyDecayAngularStrategy::execute() | "
				"Wrong parameter type!");

	size_t size = DataPointStorage::Instance().getNumberOfEvents(storage_index_);
	std::vector<std::complex<double> > values;
	values.reserve( size );
	for (unsigned int i = 0; i < size; ++i) {
		values.push_back(
				tbd_amp_->evaluate(
						storage_index_,
						i,
						paras.GetMultiUnsignedInteger(0)->GetValue(0)
				)
		);
	}
	out = std::shared_ptr<MultiComplex>(
			new MultiComplex(out->GetName(), values)
	);
	return true;
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
