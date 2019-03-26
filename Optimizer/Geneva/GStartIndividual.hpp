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

#ifndef COMPWA_OPTIMIZER_GENEVA_GSTARTINDIVIDUAL_HPP_
#define COMPWA_OPTIMIZER_GENEVA_GSTARTINDIVIDUAL_HPP_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Core/ParameterList.hpp"
#include "Estimator/Estimator.hpp"

// Global checks, defines and includes needed for all of Geneva
#include "common/GGlobalDefines.hpp"
#include "geneva/GConstrainedDoubleObject.hpp"
#include "geneva/GParameterSet.hpp"

namespace Gem {
namespace Geneva {

class GStartIndividual : public GParameterSet {
public:
  GStartIndividual(std::shared_ptr<ComPWA::Estimator::Estimator> data,
                   ComPWA::ParameterList list);

  GStartIndividual(const GStartIndividual &);
  virtual ~GStartIndividual() = default;

  bool getPar(ComPWA::ParameterList &val);

protected:
  void updatePar();

  /** @brief Loads the data of another GStartIndividual */
  void load_(const GObject *) final;

  /** @brief The actual fitness calculation takes place here. */
  double fitnessCalculation() final;

private:
  /**
   * The default constructor. Intentionally private and empty, as it is only
   * needed for serialization purposes.
   */
  GStartIndividual();

  // You can add other variables here. Do not forget to serialize them if
  // necessary int myVar;
  ComPWA::ParameterList parList;
  std::vector<std::string> parNames;
  std::shared_ptr<ComPWA::Estimator::Estimator> theData;

  /** @brief Make the class accessible to Boost.Serialization */
  friend class boost::serialization::access;

  /**
   * This function triggers serialization of this class and its
   * base classes.
   */
  template <typename Archive> void serialize(Archive &ar, const unsigned int) {
	LOG(DEBUG) << "calling serialize!!!";
    using boost::serialization::make_nvp;
    ar &BOOST_SERIALIZATION_BASE_OBJECT_NVP(GParameterSet);
    // Add other variables here like this:
    //ar &BOOST_SERIALIZATION_NVP(theData);
  }

  /** @brief Creates a deep clone of this object */
  GObject *clone_() const final;
};

} /* namespace Geneva */
} /* namespace Gem */

BOOST_CLASS_EXPORT_KEY(Gem::Geneva::GStartIndividual)

#endif
