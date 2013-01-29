#include <memory>

#include "Core/Data.hpp"
#include "Core/Estimator.hpp"
#include "Core/Amplitude.hpp"
#include "Core/Optimizer.hpp"

#include "Core/RunManager.hpp"

RunManager::RunManager(std::shared_ptr<Data> inD, std::shared_ptr<Estimator> inE,
      std::shared_ptr<Amplitude> inP, std::shared_ptr<Optimizer> inO)
      : pData_(inD), pEsti_(inE), pPhys_(inP), pOpti_(inO), valid_(false), success_(false) {
  if(inD && inE && inP && inO)
    valid_ = true;
}

RunManager::~RunManager(){
  /* nothing */
}

bool RunManager::startFit(std::vector<std::shared_ptr<PWAParameter> >& inPar){
  if( !valid_ )
    return false;

  pOpti_->exec(inPar);
  success_ = true;

  return success_;
}
