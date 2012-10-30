#ifndef PIFBASE_HPP_
#define PIFBASE_HPP_

class PIFBase
{

public:

  PIFBase()
	  {
	  }

  virtual ~PIFBase()
	{ /* nothing */	}

  virtual const double intensity(double x, double M, double T) =0;

};

#endif
