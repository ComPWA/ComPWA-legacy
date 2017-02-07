/*
 * NonResonant.cpp
 *
 *  Created on: Jan 13, 2015
 *      Author: weidenka
 */

#include "Physics/AmplitudeSum/NonResonant.hpp"

namespace ComPWA { namespace Physics { namespace AmplitudeSum {
    
    NonResonant::NonResonant(const char *name,
                             std::shared_ptr<DoubleParameter> mag,
                             std::shared_ptr<DoubleParameter> phase,
                             std::string mother, std::string particleA, std::string particleB,
                             int nCalls, normStyle nS) :
				AmpAbsDynamicalFunction(name, 0, 0, mag, phase,
                                        std::make_shared<DoubleParameter>("mass", 0.0),
                                        ComPWA::Spin(0.0), ComPWA::Spin(0.0), ComPWA::Spin(0.0),
                                        +1, 0, mother, particleA, particleB,
                                        formFactorType::noFormFactor, nCalls, nS)
    {
    }
    
    std::complex<double> NonResonant::dynamicalFunction(){
        return std::complex<double>(1,0);
    }
    
    //! Configure resonance from ptree
    void NonResonant::Configure(boost::property_tree::ptree::value_type const& v,
                                ParameterList& list)
    {
        if( v.first != "NonResonant" ) throw BadConfig("");
        
        boost::property_tree::ptree pt = v.second;
        AmpAbsDynamicalFunction::Configure(v,list);
        
        return;
    }
    
    void NonResonant::Save(boost::property_tree::ptree &pt)
    {
        boost::property_tree::ptree amp;
        AmpAbsDynamicalFunction::put(amp);
        pt.add_child("NonResonant", amp);
        return;
    }
    
    std::shared_ptr<FunctionTree> NonResonant::SetupTree(
                                                         ParameterList& sample, ParameterList& toySample,std::string suffix)
    {
        double phspVol = Kinematics::instance()->GetPhspVolume();
        
        LOG(info) << "NonResonant::setupBasicTree() | "<<_name;
        //------------Setup Tree---------------------
        std::shared_ptr<FunctionTree> newTree(new FunctionTree());
        
        int sampleSize = sample.GetMultiDouble(0)->GetNValues();
        
        //----Strategies needed
        std::shared_ptr<MultAll> mmultStrat(new MultAll(ParType::MCOMPLEX));
        std::shared_ptr<AbsSquare> msqStrat(new AbsSquare(ParType::MDOUBLE));
        std::shared_ptr<MultAll> multStrat(new MultAll(ParType::COMPLEX));
        std::shared_ptr<MultAll> multDStrat(new MultAll(ParType::DOUBLE));
        std::shared_ptr<AddAll> addStrat(new AddAll(ParType::DOUBLE));
        std::shared_ptr<Complexify> complStrat(new Complexify(ParType::COMPLEX));
        std::shared_ptr<Inverse> invStrat(new Inverse(ParType::DOUBLE));
        std::shared_ptr<SquareRoot> sqRootStrat(new SquareRoot(ParType::DOUBLE));
        
        newTree->createHead("Reso_"+_name, mmultStrat, sampleSize);
        
        newTree->createNode("PreFactor_"+_name, complStrat, "Reso_"+_name);
        newTree->createLeaf(
                            "IntensPre_"+_name, std::abs(_prefactor), "PreFactor_"+_name);
        newTree->createLeaf(
                            "PhasePre_"+_name, std::arg(_prefactor), "PreFactor_"+_name);
        
        newTree->createNode("C_"+_name, complStrat, "Reso_"+_name); //c=r*exp(phi)
        newTree->createLeaf("Intens_"+_name, _mag, "C_"+_name); //r
        newTree->createLeaf("Phase_"+_name, _phase, "C_"+_name); //phi
        
        std::shared_ptr<MultiComplex> unitVec(
                                              new MultiComplex("unit",std::vector<std::complex<double> >(
                                                                                                         sampleSize, std::complex<double>(1,0))) );
        
        newTree->createLeaf("NonRes_"+_name, unitVec, "Reso_"+_name); //nonReso
        //adding nodes and leafs for calculation of normalization
        if(_normStyle==normStyle::none){
            newTree->createLeaf("N_"+_name, 1., "Reso_"+_name);
        }else{
            newTree->createLeaf("N_"+_name, 1/std::sqrt(phspVol), "Reso_"+_name);
            //		newTree->createLeaf("N_"+_name, 1/phspVol, "Reso_"+_name);
            //		newTree->createNode("N_"+_name, sqRootStrat, "Reso_"+_name);
            //		newTree->createNode("NSq_"+_name, multDStrat, "N_"+_name);
            //		newTree->createLeaf("PhspSize_"+_name, toySampleSize, "NSq_"+_name);
            //		newTree->createLeaf("PhspVolume_"+_name, 1/phspVol, "NSq_"+_name);
            //		newTree->createNode("InvSum_"+_name, invStrat, "NSq_"+_name);
            //		newTree->createNode("Sum_"+_name, addStrat, "InvSum_"+_name);
            //		newTree->createNode("AbsVal_"+_name, msqStrat, "Sum_"+_name);
            //		std::shared_ptr<MultiComplex> unitVec2(
            //				new MultiComplex("unit",std::vector<std::complex<double> >(
            //						toySampleSize, std::complex<double>(1,0))) );
            //		newTree->createLeaf("NormNonRes_"+_name, unitVec2, "AbsVal_"+_name);
        }
        return newTree;
    }
} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */
