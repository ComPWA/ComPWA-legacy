//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
#include <memory>

#include "DataReader/Data.hpp"
#include "Estimator/Estimator.hpp"
#include "Physics/Amplitude.hpp"
#include "Optimizer/Optimizer.hpp"
#include "Core/ParameterList.hpp"

#include "Core/RunManager.hpp"

RunManager::RunManager(std::shared_ptr<Data> inD, std::shared_ptr<ControlParameter> inE,
      std::shared_ptr<Amplitude> inP, std::shared_ptr<Optimizer> inO)
      : pData_(inD), pEsti_(inE), pPhys_(inP), pOpti_(inO), valid_(false), success_(false) {
  if(inD && inE && inP && inO)
    valid_ = true;
}

RunManager::~RunManager(){
  /* nothing */
}

bool RunManager::startFit(ParameterList& inPar){
  if( !valid_ )
    return false;

  pOpti_->exec(inPar);
  success_ = true;

  return success_;
}
