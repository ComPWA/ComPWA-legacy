#include <memory>

#include "DIFBase.hpp"
#include "EIFBase.hpp"
#include "PIFBase.hpp"
#include "OIFBase.hpp"

#include "RunManager.hpp"

RunManager::RunManager(std::shared_ptr<DIFBase> inD, std::shared_ptr<EIFBase> inE,
      std::shared_ptr<PIFBase> inP, std::shared_ptr<OIFBase> inO)
      : pData_(inD), pEsti_(inE), pPhys_(inP), pOpti_(inO), valid_(false), success_(false) {
  if(inD && inE && inP && inO)
    valid_ = true;
}

RunManager::~RunManager(){
  /* nothing */
}

bool RunManager::startFit(std::vector<PWAParameter<double> >& inPar){
  if( !valid_ )
    return false;

  pOpti_->exec(inPar);
  success_ = true;

  return success_;
}
