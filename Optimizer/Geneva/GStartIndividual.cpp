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

namespace Gem {
namespace Geneva {

GStartIndividual::GStartIndividual() : GParameterSet(), theData() {}

GStartIndividual::GStartIndividual(
    std::shared_ptr<ComPWA::Estimator::Estimator> data,
    ComPWA::ParameterList list)
    : GParameterSet(), parList(list), theData(data) {
  for (std::size_t i = 0; i < parList.doubleParameters().size(); i++) {
    auto p = parList.doubleParameter(i);
    if (p->isFixed())
      continue;
    double val = p->value();
    double min = GConstrainedValueLimitT<double>::lowest();
    double max = GConstrainedValueLimitT<double>::highest();
    double err = val;
    if (p->hasError())
      err = p->error().first;
    if (p->hasBounds()) {
      min = p->bounds().first;
      max = p->bounds().second;
    }
    std::shared_ptr<GConstrainedDoubleObject> gbd_ptr(
        new GConstrainedDoubleObject(val, min, max));

    // TODO: I don't understand what these values mean. I simply chose some. Fix
    // this!
    std::shared_ptr<GDoubleGaussAdaptor> gdga_ptr(
        new GDoubleGaussAdaptor(0.5, 0.5, 0., 0.95));
    gdga_ptr->setAdaptionThreshold(
        1); // Adaption parameters are modified after each adaption
    gdga_ptr->setAdaptionProbability(
        0.05); // The likelihood for a parameter to be adapted

    // Register the adaptor with GConstrainedDoubleObject objects
    gbd_ptr->addAdaptor(gdga_ptr);

    this->push_back(gbd_ptr);
    parNames.push_back(p->name());
  }
  LOG(INFO) << "GStartIndividual::GStartIndividual() | " << parNames.size()
            << " Parameters were added for minimization!";
}

GStartIndividual::GStartIndividual(const GStartIndividual &cp)
    : GParameterSet(cp), parNames(cp.parNames), theData(cp.theData) {
  // TODO: the deep copy is necessary to achieve thread safety, however I
  parList.DeepCopy(cp.parList);
}

bool GStartIndividual::getPar(ComPWA::ParameterList &val) {
  updatePar();
  val = ComPWA::ParameterList(parList);
  return true;
}

void GStartIndividual::load_(const GObject *cp) {
  const GStartIndividual *p_load =
      Gem::Common::g_convert_and_compare<GObject, GStartIndividual>(cp, this);
  // Load our parent's data
  GParameterSet::load_(cp);

  // TODO: I'm not sure if you are supposed to load the Estimator or Parameters
  // here
}

GObject *GStartIndividual::clone_() const {
  return new GStartIndividual(*this);
}

double GStartIndividual::fitnessCalculation() {
  updatePar();
  return theData->evaluate();
}

void GStartIndividual::updatePar() {
  // Somehow you need to set this->end() first. Setting the iterator to
  // this->begin() right away, will not work...
  GStartIndividual::conversion_iterator<GConstrainedDoubleObject> it(
      this->end());
  it = this->begin();
  for (unsigned int i = 0; i < parList.doubleParameters().size(); ++i) {
    auto p = parList.doubleParameter(i);
    if (p->isFixed())
      continue;
    p->setValue((*it)->value());
    ++it;
  }
}
} /* namespace Geneva */
} /* namespace Gem */
