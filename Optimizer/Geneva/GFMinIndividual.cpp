// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

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

#include "GFMinIndividual.hpp"

#include "Core/Logging.hpp"

#include "geneva/GConstrainedDoubleObject.hpp"
#include "geneva/GDoubleGaussAdaptor.hpp"

BOOST_CLASS_EXPORT_IMPLEMENT(Gem::Geneva::GFMinIndividual)

namespace Gem {
namespace Geneva {

GFMinIndividual::GFMinIndividual(
    ComPWA::Estimator::Estimator<double> &estimator,
    ComPWA::FitParameterList parlist)
    : Estimator(&estimator) {
  auto ActualParameters = estimator.getParameters();

  size_t varcounter(0);
  // it is assumed that the Estimator Parameters and FitParameterList are
  // synchronized (this was validated previously)
  for (auto const &p : parlist) {
    if (!p.IsFixed) {
      FreeParameterIndices.push_back(varcounter);
      double val = p.Value;
      double min = GConstrainedValueLimitT<double>::lowest();
      double max = GConstrainedValueLimitT<double>::highest();

      if (p.HasBounds) {
        min = p.Bounds.first;
        max = p.Bounds.second;
      }
      std::shared_ptr<GConstrainedDoubleObject> gbd_ptr(
          new GConstrainedDoubleObject(val, min, max));

      std::shared_ptr<GDoubleGaussAdaptor> gdga_ptr(
          new GDoubleGaussAdaptor(GFI_DEF_SIGMA, GFI_DEF_SIGMASIGMA,
                                  GFI_DEF_MINSIGMA, GFI_DEF_MAXSIGMA));
      gdga_ptr->setAdaptionThreshold(
          1); // Adaption parameters are modified after each adaption
      gdga_ptr->setAdaptionProbability(
          GFI_DEF_ADPROB); // The likelihood for a parameter to be adapted

      // Register the adaptor with GConstrainedDoubleObject objects
      gbd_ptr->addAdaptor(gdga_ptr);
      this->push_back(gbd_ptr);
    }
    AllParameters.push_back(p.Value);
    ++varcounter;
  }
  LOG(INFO) << "GStartIndividual::GStartIndividual() | "
            << FreeParameterIndices.size()
            << " Parameters were added for minimization!";
}

GFMinIndividual::GFMinIndividual(const GFMinIndividual &cp)
    : GParameterSet(cp), Estimator(cp.Estimator),
      AllParameters(cp.AllParameters),
      FreeParameterIndices(cp.FreeParameterIndices) {}

const GFMinIndividual &GFMinIndividual::operator=(const GFMinIndividual &cp) {
  GFMinIndividual::load_(&cp);
  return *this;
}

void GFMinIndividual::load_(const GObject *cp) {
  // Check that we are dealing with a GFMinIndividual reference independent of
  // this object and convert the pointer
  const GFMinIndividual *p_load =
      Gem::Common::g_convert_and_compare<GObject, GFMinIndividual>(cp, this);

  // Load our parent class'es data ...
  GParameterSet::load_(cp);

  // ... and then our local data
  Estimator = p_load->Estimator;
  AllParameters = p_load->AllParameters;
  FreeParameterIndices = p_load->FreeParameterIndices;
}

GObject *GFMinIndividual::clone_() const { return new GFMinIndividual(*this); }

double GFMinIndividual::fitnessCalculation() {
  // Retrieve the parameters
  std::vector<double> NewFreeParameters;
  this->streamline(NewFreeParameters);

  auto NewParameters = AllParameters;
  for (unsigned int i = 0; i < NewFreeParameters.size(); ++i) {
    NewParameters[FreeParameterIndices[i]] = NewFreeParameters[i];
  }

  Estimator->updateParametersFrom(NewParameters);
  return Estimator->evaluate();
}

} // namespace Geneva
} // namespace Gem
