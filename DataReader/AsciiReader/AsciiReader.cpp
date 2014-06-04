//-------------------------------------------------------------------------------
// Copyright (c) 2013 Florian Feldbauer.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Florian Feldbauer - initial API and implementation
//-------------------------------------------------------------------------------
//_____ I N C L U D E S _______________________________________________________

// ANSI C headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>

// 3rd party headers
//#include "external/qft++/include/Tensor.h"

// local headers
#include "Core/Exceptions.hpp"
#include "DataReader/AsciiReader/AsciiReader.h"

//_____ D E F I N I T I O N S __________________________________________________

//_____ G L O B A L S __________________________________________________________

//_____ L O C A L S ____________________________________________________________

//_____ F U N C T I O N S ______________________________________________________

// Constructors and destructors
AsciiReader::AsciiReader( const std::string inConfigFile, const int particles  )
  : fmaxBins_(0)
{

  std::ifstream currentStream;
  currentStream.open ( inConfigFile.c_str() );

  if ( !currentStream ) throw BadConfig( "Can not open " + inConfigFile );

  while( !currentStream.eof() ) {
    double e, px, py, pz;
    Event newEvent;

    for ( int parts = 0; parts < particles; parts++ ) {
      currentStream >> px >> py >> pz >> e;
      newEvent.addParticle( Particle( px, py, pz, e ) );
    }

    if (!currentStream.fail()) {
      EvtList_.push_back( newEvent );
      // for ( parts = 0; parts < linesToSkip; parts++ )
      //   currentStream >> px >> py >> pz >> e;
    }
  }
  currentStream.close();
}

AsciiReader::~AsciiReader() {
  EvtList_.clear();
}

allMasses AsciiReader::getMasses(){
  if(!EvtList_.size()) return allMasses();
  unsigned int nParts = EvtList_.at(0).getNParticles();
  BOOST_LOG_TRIVIAL(debug)<<"RootReader::getMasses() #particles: "<<nParts;

  //determine invMass combinations
  unsigned int nMasses=0;
  std::vector<std::pair<unsigned int, unsigned int> > ids;
  for(unsigned int i=0; i<nParts; i++)
    for(unsigned int j=i+1; j<nParts; j++){
      nMasses++;
      ids.push_back(std::make_pair(i+1,j+1));
    }
  BOOST_LOG_TRIVIAL(debug)<<"AsciiReader::getMasses() #invMasses: "<<nMasses;

  allMasses result(nMasses, EvtList_.size(), ids);
  //calc and store inv masses
  for(unsigned int evt=0; evt<EvtList_.size(); evt++){
    Event tmp = EvtList_.at(evt);

    // Check number of particle in TClonesrray
    if( nParts != tmp.getNParticles() ){
      result.nEvents--;
      continue;
    }

    for(unsigned int pa=0; pa<nParts; pa++){
      for(unsigned int pb=pa+1; pb<nParts; pb++){
        const Particle &inA(tmp.getParticle(pa));
        const Particle &inB(tmp.getParticle(pb));
        double mymass_sq = inA.invariantMass(inB);

        (result.masses_sq.at(std::make_pair(pa+1,pb+1))).at(evt) = mymass_sq;

        //tmp.addParticle(Particle(inN.X(), inN.Y(), inN.Z(), inN.E(),partN->GetPdgCode()));
        //tmp.setWeight(feventWeight); //Todo: weight? what weight? lets wait...


      }//particle loop B
    }//particle loop A

  }//event loop

  return result;
}
Event& AsciiReader::getEvent( const int index ) {
  if ( EvtList_.size() <= (unsigned int)index )
    throw BadIndex("Index exceeds max number of events");

  return EvtList_.at( index );
}

const int AsciiReader::getBin( const int i, double& m12, double& weight ) {
  return 0;
}

