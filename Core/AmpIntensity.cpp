/*
 * AmpIntensity.cpp
 *
 *  Created on: Mar 16, 2016
 *      Author: weidenka
 */

#include "Core/AmpIntensity.hpp"

namespace ComPWA {

void AmpIntensity::UpdateParameters(ParameterList &par) {
  std::shared_ptr<DoubleParameter> pOld, pNew;

  /* First we check if new parameter list contains at least one matching
   * parameter. Otherwise we skip! */
  int commonPar = 0;
  for (unsigned int i = 0; i < params.GetNDouble(); i++) {
    try {
      pNew = par.GetDoubleParameter(params.GetDoubleParameter(i)->GetName());
    } catch (std::exception &ex) {
      continue;
    }
    commonPar++;
  }
  if (commonPar == 0)
    return;

  /* If we have at least one matching parameter, we require that all
   * parameters are contained in the new list */
  for (unsigned int i = 0; i < params.GetNDouble(); i++) {
    pOld = params.GetDoubleParameter(i);
    try {
      pNew = par.GetDoubleParameter(pOld->GetName());
    } catch (std::exception &ex) {
      LOG(error) << "AmpSumIntensity::setParameterList() | "
                    " Can not find parameter! "
                 << ex.what();
      throw;
    }
    // Update parameter
    pOld->UpdateParameter(pNew);
  }

  return;
}

void AmpIntensity::FillParameterList(ParameterList &outPar) const {
  // Parameters are only added if they do not exist yet
  outPar.Append(params);
  outPar.RemoveDuplicates();
}
  
} /* namespace ComPWA */
