#ifndef _PIFBW_HPP
#define _PIFBW_HPP

#include <vector>
#include <memory>
#include "PIFBase.hpp"


class PIFBW : public PIFBase {

public:
  /// Default Constructor (0x0)
  PIFBW();

  virtual const double intensity(double x, double M, double T);

  /** Destructor */
  virtual ~PIFBW();

 protected:

 private:
  const double BreitWigner(double x, double M, double T);

};

#endif /* _PIFBW_HPP */
