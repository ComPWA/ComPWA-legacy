//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - correct nominator, using dataPoint for data handling
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#ifndef AMP_REL_BREIT_WIGNER_RES
#define AMP_REL_BREIT_WIGNER_RES

#include <vector>

//#include "TObject.h"
//#include "TString.h"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"
#include "Physics/AmplitudeSum/AmpKinematics.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"
#include "Core/Functions.hpp"
#include "Core/Exceptions.hpp"

//#include "Physics/AmplitudeSum/AmpWigner.hpp"

class BreitWignerStrategy : public Strategy {
public:
  BreitWignerStrategy(const std::string resonanceName):name(resonanceName){
    //name = +resonanceName;
  }

  virtual const std::string to_str() const {
    return ("relativistic BreitWigner of "+name);
  }

  virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter> out) {
    out = std::shared_ptr<AbsParameter>();

    double Gamma0, GammaV, m0, d, norm, BLWeiss2, qTerm;
    unsigned int spin, subSys;
    try{
      m0 = double(paras.GetParameterValue("m0_"+name));
    }catch(BadParameter& e){
      m0 = double(paras.GetParameterValue("ParOfNode_m0_"+name));
    }
    spin = (unsigned int)(paras.GetParameterValue("ParOfNode_spin_"+name));
    d = double(paras.GetParameterValue("ParOfNode_d_"+name));
    norm = double(paras.GetParameterValue("ParOfNode_norm_"+name));
    subSys = double(paras.GetParameterValue("ParOfNode_subSysFlag_"+name));

    //m  = double(paras.GetParameterValue("mym")); //TODO: MultiDIM!!!
   // ma = double(paras.GetParameterValue("ma"));
   // mb = double(paras.GetParameterValue("mb"));

    try{
      Gamma0 = double(paras.GetParameterValue("resWidth_"+name));
    }catch(BadParameter& e){
      Gamma0 = double(paras.GetParameterValue("ParOfNode_resWidth_"+name));
    }
    /*GammaV = Gamma0 * (m0 / m) * pow(q(ma,mb,m) / q0(ma,mb,m0), 2.*spin + 1.)  * BLprime2(ma,mb,m0,m,d,spin);

    std::complex<double> denom(m0*m0 - m*m, -m0 * GammaV);
    std::complex<double> res(m0 * Gamma0);
    res = res / denom;*/

    //MultiDim output, must have multidim Paras in input
    if(out->type() == ParType::MCOMPLEX){
      if(paras.GetNMultiDouble()){
        unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();

        std::vector<std::complex<double> > results(nElements, std::complex<double>(0.));
        std::shared_ptr<MultiDouble> mp, map, mbp;
        switch(subSys){
          case 3:{ //reso in sys of particles 1&2
            mp  = (paras.GetMultiDouble("m12"));
            map  = (paras.GetMultiDouble("m23"));
            mbp  = (paras.GetMultiDouble("m13"));
            break;
          }
          case 4:{ //reso in sys of particles 1&3
            mp  = (paras.GetMultiDouble("m13"));
            map  = (paras.GetMultiDouble("m12"));
            mbp  = (paras.GetMultiDouble("m23"));
            break;
          }
          case 5:{ //reso in sys of particles 2&3
            mp  = (paras.GetMultiDouble("m23"));
            map  = (paras.GetMultiDouble("m13"));
            mbp  = (paras.GetMultiDouble("m12"));
            break;
          }
        }

        //calc BW for each point
        for(unsigned int ele=0; ele<nElements; ele++){
          BLWeiss2 = BLprime2(map->GetValue(ele),mbp->GetValue(ele),m0,mp->GetValue(ele),d,spin);
          qTerm = pow(q(map->GetValue(ele),mbp->GetValue(ele),mp->GetValue(ele)) / q0(map->GetValue(ele),mbp->GetValue(ele),m0), 2.*spin + 1.);
          //double Gamma0 = _resWidth.GetValue();
          GammaV = Gamma0 * qTerm * (m0 / mp->GetValue(ele)) * BLWeiss2;
          std::complex<double> denom(m0*m0 - mp->GetValue(ele)*mp->GetValue(ele), -m0 * GammaV);

          results[ele] = (std::complex<double>(norm) / denom); //Laura++ (old) definition*/
        }

        out = std::shared_ptr<AbsParameter>(new MultiComplex("relBW of "+name,results));
        return true;
      }else{ //end multidim para treatment
        //Todo: exception wrong input
        return false;
      }
    }//end multicomplex output


    //Only StandardDim Paras in input
    //  double spinTerm = evaluateWignerD(); //spinTerm =1;
    double m, ma, mb;
    switch(subSys){
      case 3:{ //reso in sys of particles 1&2
        m  = double(paras.GetParameterValue("m12"));
        ma  = double(paras.GetParameterValue("m23"));
        mb  = double(paras.GetParameterValue("m13"));
        break;
      }
      case 4:{ //reso in sys of particles 1&3
        m  = double(paras.GetParameterValue("m13"));
        ma  = double(paras.GetParameterValue("m12"));
        mb  = double(paras.GetParameterValue("m23"));
        break;
      }
      case 5:{ //reso in sys of particles 2&3
        m  = double(paras.GetParameterValue("m23"));
        ma  = double(paras.GetParameterValue("m13"));
        mb  = double(paras.GetParameterValue("m12"));
        break;
      }
    }
    BLWeiss2 = BLprime2(ma,mb,m0,m,d,spin);
    qTerm = pow(q(ma,mb,m) / q0(ma,mb,m0), 2.*spin + 1.);
    //double Gamma0 = _resWidth.GetValue();
    GammaV = Gamma0 * qTerm * (m0 / m) * BLWeiss2;
    std::complex<double> denom(m0*m0 - m*m, -m0 * GammaV);

    std::complex<double> result = std::complex<double>(norm) / denom; //Laura++ (old) definition*/

    //std::complex<double> result (res.re(),res.im());
    out = std::shared_ptr<AbsParameter>(new ComplexParameter("relBW of "+name, result));
    return true;
  }

protected:
  std::string name;

  double q0(const double& ma, const double& mb, const double& m0) const {
    double mapb = ma + mb;
    double mamb = ma - mb;

    return sqrt ( (m0*m0 - mapb*mapb) * (m0*m0 - mamb*mamb) ) / (2. * m0 );
  }

  double q(const double& ma, const double& mb, const double& x) const {
    double mapb = ma + mb;
    double mamb = ma - mb;

    return sqrt ( (x*x - mapb*mapb) * (x*x - mamb*mamb) ) / (2. * x );
  }


  // compute part of the Blatt-Weisskopf barrier factor
  //   BLprime = sqrt (F(q0)/F(q))
  double F(const double& p, const double& d, unsigned int& spin) const {
    double retVal = 1;

    if (spin == 0)
      retVal = 1;
    else if (spin == 1)
      retVal = 1 + p*p * d*d;
    else if (spin == 2) {
      double z = p*p * d*d;
      retVal = (z-3.)*(z-3.) + 9*z;
    }
    return retVal;
  }


  // compute square of Blatt-Weisskopf barrier factor
  double BLprime2(const double& ma, const double& mb, const double& m0, const double& x, const double& d, unsigned int& spin) const {
    //  cout << q0() << " " << q() << "\t" << F(q0()) << " " << F(q()) << endl;
    return F(q0(ma, mb, m0),d,spin) / F(q(ma, mb, x),d,spin);
  }

};

class AmpRelBreitWignerRes : public AmpAbsDynamicalFunction, public AmpKinematics {
public:

	AmpRelBreitWignerRes(const char *name,
			DoubleParameter& _resMass, DoubleParameter& _resWidth,
			double& _radius,
			int _subsys,
			int resSpin, int m, int n
	) ;
	AmpRelBreitWignerRes(const AmpRelBreitWignerRes&, const char*);
	AmpRelBreitWignerRes(const AmpRelBreitWignerRes&);

	virtual ~AmpRelBreitWignerRes();

	//  double operator() (double *x, size_t dim, void*);

	virtual void initialise();
	virtual std::complex<double> evaluate(dataPoint& point) ;
	virtual std::complex<double> evaluateAmp(dataPoint& point);
	virtual double evaluateWignerD(dataPoint& point) { return _wignerD.evaluate(point); };
	//  virtual double eval(double x[],size_t dim, void *param) const;//used for MC integration
	//  double (*eval2)(double x[],size_t dim, void *param);//used for MC integration

	void setDecayMasses(double, double, double, double);
	//  double getMaximum() const{return 1;};
	double getSpin() {return _spin;}; //needs to be declared in AmpAbsDynamicalFunction
	inline virtual bool isSubSys(const unsigned int subSys) const{ return (subSys==_subSys); };

protected:
	DoubleParameter _resWidth;
	AmpWigner2 _wignerD;
	bool foundMasses;
	unsigned int id23, id13;
//	AmpWigner _wignerD;

private:

};

#endif
