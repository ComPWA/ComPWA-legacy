//! Physics Module with simple 1D Breit-Wigner.
/*! \class PIFBW
 * @file PIFBW.hpp
 * This class provides a simple Breit-Wigner calculation with given parameters.
 * It fulfills the Amplitude interface to be compatible with other ComPWA modules.
*/

#ifndef _PIFBW_HPP
#define _PIFBW_HPP

#include <vector>
#include <memory>
#include "Amplitude.hpp"


class PIFBW : public Amplitude {

public:
  /// Default Constructor (0x0)
  PIFBW(const double min, const double max);

  //For normalization
  virtual const double integral(std::vector<std::shared_ptr<PWAParameter> >& par);
  virtual const double drawInt(double *x, double *p); //For easy usage in a root TF1
  virtual const double intensity(double x, double M, double T);
  virtual const double intensity(std::vector<double> x, std::vector<std::shared_ptr<PWAParameter> >& par);
  virtual const bool fillStartParVec(std::vector<std::shared_ptr<PWAParameter> >& outPar);

  /** Destructor */
  virtual ~PIFBW();

protected:
  double min_;
  double max_;

 private:
  const double BreitWigner(double x, double M, double T);

};

#endif /* _PIFBW_HPP */
