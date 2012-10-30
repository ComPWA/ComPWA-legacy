#ifndef MIBASE_HPP_
#define MIBASE_HPP_

class OIFBase
{

public:

  OIFBase()
	  {
	  }

  virtual ~OIFBase()
	{ /* nothing */	}

  virtual const double exec(unsigned int Npar, double* par,  double* min, double* max, double* err) =0;
 
};

#endif
