//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#ifndef AMP_GAUS_RES
#define AMP_GAUS_RES

#include <vector>

#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"

namespace ComPWA { namespace Physics { namespace AmplitudeSum {
    
    class AmpGausRes : public AmpAbsDynamicalFunction  {
    public:
        
        AmpGausRes();
        AmpGausRes(const char *name,
                   unsigned int varIdA,
                   std::shared_ptr<DoubleParameter> mag,
                   std::shared_ptr<DoubleParameter> phase,
                   std::shared_ptr<DoubleParameter> mass,
                   std::shared_ptr<DoubleParameter> width,
                   std::string mother, std::string particleA, std::string particleB,
                   int nCalls=30000, normStyle nS=normStyle::one) ;
        
        //! Clone function
        virtual AmpGausRes* Clone(std::string newName="") const{
            auto tmp = (new AmpGausRes(*this));
            if(newName != "")
                tmp->SetName(newName);
            return tmp;
        }
        
        ~AmpGausRes();
        
        //! Calculation integral |dynamical amplitude|^2
        virtual double GetIntegral() { return integral(); }
        
        //! Get resonance width
        double GetWidth() const { return _width->GetValue(); }
        
        virtual void Save(boost::property_tree::ptree&) { };
        
        virtual std::complex<double> Evaluate(dataPoint& point);
        virtual std::complex<double> EvaluateAmp(dataPoint& point);
        virtual double evaluateWignerD(dataPoint& point) const { return 1; };
        
        inline virtual bool isSubSys(const unsigned int subSys) const
        {
            return (subSys==_subSys);
        };
        
        double GetSpin() const { return 0; };
        
        virtual std::shared_ptr<FunctionTree>
        SetupTree(ParameterList& sample,
                  ParameterList& toySample,
                  std::string suffix) {
            return std::shared_ptr<FunctionTree>();
        };
        
    protected:
        std::shared_ptr<DoubleParameter> _width;
        
    };
    
} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */

#endif
