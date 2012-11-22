//! Run-Manager for a simple fit.
/*! \class RunManager
 * @file RunManager.hpp
 * This class provides a RunManager for simple fits. To use it, you create
 * all modules you want to use and provide them to the RunManger. It checks
 * for compatibility and if set up correctly it starts the fitting procedure.
*/

#ifndef _RUNMANAGER_HPP_
#define _RUNMANAGER_HPP_

#include <vector>
#include <memory>

#include "DIFBase.hpp"
#include "EIFBase.hpp"
#include "PIFBase.hpp"
#include "OIFBase.hpp"
#include "PWAParameter.hpp"

class RunManager
{

public:

  RunManager(std::shared_ptr<DIFBase>, std::shared_ptr<EIFBase>,
      std::shared_ptr<PIFBase>, std::shared_ptr<OIFBase>);

  virtual ~RunManager();

  virtual bool startFit(std::vector<std::shared_ptr<PWAParameter> >& );

protected:
  std::shared_ptr<DIFBase> pData_; /*!< Pointer to Data-Module */
  std::shared_ptr<EIFBase> pEsti_; /*!< Pointer to Estimator-Module */
  std::shared_ptr<PIFBase> pPhys_; /*!< Pointer to Physics-Module */
  std::shared_ptr<OIFBase> pOpti_; /*!< Pointer to Optimizer-Module */
  //TODO: log
  bool valid_; /*!< setup a valid configuration? */
  bool success_; /*!< fitting ended successfully? */

};

#endif
