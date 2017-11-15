//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding functionality to generate set of
//events
//-------------------------------------------------------------------------------
//! Run-Manager for a simple fit.
/*! \class PythonFit
 * @file PythonFit.hpp
 * This class provides a Manager for simple fits. It creates a set of modules
 * for an unbinned likelihood fit of a three-body final state.
 */

#ifndef _PYTHONFIT_HPP_
#define _PYTHONFIT_HPP_

#include <vector>
#include <memory>
#include <string>

#include <boost/program_options.hpp>
#include <boost/serialization/export.hpp>

#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TPython.h"

// Core header files go here
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/TableFormater.hpp"
#include "Core/AbsParameter.hpp"
#include "Core/Logging.hpp"

// ComPWA header files go here
#include "DataReader/RootReader/RootReader.hpp"
#include "DataReader/RootReader/RootEfficiency.hpp"
#include "DataReader/CorrectionTable.hpp"
#include "DataReader/DataCorrection.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Optimizer/Minuit2/MinuitResult.hpp"
//#include "Optimizer/Geneva/GenevaIF.hpp"
//#include "Optimizer/Geneva/GenevaResult.hpp"
#include "Physics/HelicityFormalism.hpp"

#include "Tools/RunManager.hpp"
#include "Tools/RootGenerator.hpp"
#include "Tools/FitFractions.hpp"

// We define an intensity model using a raw string literal. Currently, this is
// just a toy model without any physical meaning.
// (comments within the string are ignored!). This is convenient since we
// do not have to configure the build system to copy input files somewhere.
// In practise you may want to use a normal XML input file instead.
std::string amplitudeModel = R"####(
<IncoherentIntensity Name="jpsiGammaPiPi_inc">
	<CoherentIntensity Name="jpsiGammaPiPi">
  	<Amplitude Name="f2(1270)">
			<Parameter Type="Magnitude"	Name="Magnitude_f2">
				<Value>1.0</Value>
        <Min>-1.0</Min>
        <Max>2.0</Max>
				<Fix>false</Fix>
			</Parameter>
			<Parameter Type="Phase" Name="Phase_f2">
				<Value>0.0</Value>
        <Min>-100</Min>
        <Max>100</Max>
				<Fix>false</Fix>
			</Parameter>
			<Resonance Name="f2ToPiPi">
				<DecayParticle Name="f2(1270)" Helicity="0"/>
				<SubSystem>
					<RecoilSystem FinalState="0" />
					<DecayProducts>
						<Particle Name="pi0" FinalState="1"  Helicity="0"/>
						<Particle Name="pi0" FinalState="2"  Helicity="0"/>
					</DecayProducts>
				</SubSystem>
			</Resonance>
		</Amplitude>
		<Amplitude Name="myAmp">
			<Parameter Type="Magnitude"	Name="Magnitude_my">
				<Value>1.0</Value>
        <Min>-1.0</Min>
        <Max>2.0</Max>
				<Fix>true</Fix>
			</Parameter>
			<Parameter Type="Phase" Name="Phase_my`">
				<Value>0.0</Value>
        <Min>-100</Min>
        <Max>100</Max>
				<Fix>true</Fix>
			</Parameter>
			<Resonance Name="MyResToPiPi">
				<DecayParticle Name="myRes" Helicity="0"/>
				<SubSystem>
					<RecoilSystem FinalState="0" />
					<DecayProducts>
						<Particle Name="pi0" FinalState="1"  Helicity="0"/>
						<Particle Name="pi0" FinalState="2"  Helicity="0"/>
					</DecayProducts>
				</SubSystem>
			</Resonance>
		</Amplitude>
	</CoherentIntensity>
</IncoherentIntensity>
)####";

std::string myParticles = R"####(
<ParticleList>
	<Particle Name="f2(1270)">
		<Pid>225</Pid>
		<Parameter Type="Mass" Name="Mass_f2(1270)">
			<Value>1.2755</Value>
			<Error>8.0E-04</Error>
      <Min>0.1</Min>
      <Max>2.0</Max>
      <Fix>false</Fix>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="2"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="+1"/>
		<QuantumNumber Class="Int" Type="Cparity" Value="+1"/>
		<DecayInfo Type="relativisticBreitWigner">
			<FormFactor Type="0" />
			<Parameter Type="Width" Name="Width_f2(1270)">
				<Value>0.1867</Value>
			</Parameter>
			<Parameter Type="MesonRadius" Name="Radius_rho">
				<Value>2.5</Value>
				<Fix>true</Fix>
			</Parameter>
		</DecayInfo>
	</Particle>
	<Particle Name="myRes">
		<Pid>999999</Pid>
		<Parameter Type="Mass" Name="Mass_myRes">
			<Value>2.0</Value>
			<Error>8.0E-04</Error>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="1"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="+1"/>
		<QuantumNumber Class="Int" Type="Cparity" Value="+1"/>
		<DecayInfo Type="relativisticBreitWigner">
			<FormFactor Type="0" />
			<Parameter Type="Width" Name="Width_myRes">
				<Value>1.0</Value>
        <Min>0.1</Min>
        <Max>1.0</Max>
        <Fix>false</Fix>
			</Parameter>
			<Parameter Type="MesonRadius" Name="Radius_myRes">
				<Value>2.5</Value>
				<Fix>true</Fix>
			</Parameter>
		</DecayInfo>
	</Particle>
</ParticleList>
)####";

/*template<typename T>
PyObject* convertToPy(const T& cxxObj) {
	T* newCxxObj = new T(cxxObj);
	return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName(), true);
}

//template<typename T>
//PyObject* convertToPy(const T* cxxObj) {
//	return TPython::ObjectProxy_FromVoidPtr(&cxxObj, cxxObj->ClassName());
//};

template<typename T>
T convertFromPy(PyObject* pyObj) {
	TObject* TObj = (TObject*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
	T cxxObj = dynamic_cast<T>(TObj);
	return cxxObj;
}

PyObject* convertTreeToPy(TTree* cxxTree) {
	TTree* newCxxTree = cxxTree->CloneTree();
	return TPython::ObjectProxy_FromVoidPtr(newCxxTree, newCxxTree->ClassName(), true);
};*/

///
/// Simple Dalitz plot fit of the channel J/psi -> gamma pi0 pi0
///
/// The basic procedure is the following:
/// 1) Create Kinematics object
/// 2) Generate a large phase space sample
/// 3) Create intensity from pre-defined model
/// 4) Generate a data sample given intensity and kinematics
/// 5) Fit the model to the data
/// 6) Plot data sample and intensity
///
class PythonFit {
public:
  PythonFit();

  virtual ~PythonFit();

  virtual int StartFit();

protected:
  int argc;
  char ** argv;


};

#endif
