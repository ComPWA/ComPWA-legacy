/*
 * Copyright (C) Gemfony scientific UG (haftungsbeschraenkt)
 *
 * See the AUTHORS file in the top-level directory for a list of authors.
 *
 * Contact: contact [at] gemfony (dot) eu
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
 * http://www.gemfony.eu .
 */

#ifndef COMPWA_OPTIMIZER_GENEVA_GFMININDIVIDUAL_HPP_
#define COMPWA_OPTIMIZER_GENEVA_GFMININDIVIDUAL_HPP_

#include "Core/FitParameter.hpp"
#include "Estimator/Estimator.hpp"

// Global checks, defines and includes needed for all of Geneva
#include "common/GGlobalDefines.hpp"

#include "geneva/GParameterSet.hpp"

namespace Gem {
namespace Geneva {

/**
 * This individual searches for a minimum in the given Estimator
 * each capable of processing their input in multiple dimensions.
 */
class GFMinIndividual : public GParameterSet {
public:
  GFMinIndividual() = default;
  GFMinIndividual(ComPWA::Estimator::Estimator<double> &estimator,
                  ComPWA::FitParameterList parlist);
  GFMinIndividual(const GFMinIndividual &);
  virtual ~GFMinIndividual() = default;

  const GFMinIndividual &operator=(const GFMinIndividual &);

protected:
  /** @brief Loads the data of another GFMinIndividual */
  void load_(const GObject *) final;

  /** @brief The actual value calculation takes place here */
  double fitnessCalculation() final;

private:
  friend class boost::serialization::access;

  template <class Archive> void serialize(Archive &ar, const unsigned int) {
    ar &BOOST_SERIALIZATION_BASE_OBJECT_NVP(GParameterSet);
  }

  /** @brief Creates a deep clone of this object */
  GObject *clone_() const final;

  // default settings for the adaption
  // TODO: overwrite these values from the user
  const double GFI_DEF_ADPROB = 0.1;
  const double GFI_DEF_SIGMA = 0.25;
  const double GFI_DEF_SIGMASIGMA = 0.8;
  const double GFI_DEF_MINSIGMA = 0.001;
  const double GFI_DEF_MAXSIGMA = 1;

  // TODO: this should not be a pointer, but rather a smart pointer
  // (references cannot be used here, because of the default constructor
  // requirement)
  ComPWA::Estimator::Estimator<double> *Estimator = nullptr;
  std::vector<double> AllParameters;
  std::vector<unsigned int> FreeParameterIndices;
};

} /* namespace Geneva */
} /* namespace Gem */

BOOST_CLASS_EXPORT_KEY(Gem::Geneva::GFMinIndividual)

#endif
