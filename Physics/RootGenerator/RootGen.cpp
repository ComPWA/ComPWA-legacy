// Root header files go here
#include "TMath.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TGenPhaseSpace.h"
#include "TFile.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TRandom3.h"


RootGen::RootGen(double inM, std::vector<double> fiM, std::vector<int> pids)
: initialM_(inM), finalM_(fiM), pids_(pids){
  if(finalM_.size()!=pids_.size()) inM=0; //todo exception;
}

RootGen::~RootGen(){

}


void RootGen::generate(unsigned int nEvents, std::shared_pointer<Amplitude> model, std::string filename){
  unsigned int evt=0;
  TRandom3 rando;

  //Output File setup
  TFile output(filename,"recreate");
  output.SetCompressionLevel(1); //try level 2 also
  TTree fTree ("data","RootGen");
  TClonesArray *fEvt = new TClonesArray("TParticle");
  TClonesArray &ar = *fEvt;
  fTree.Branch("Particles",&fEvt);

  //initial energy:
  TLorentzVector W(0.0, 0.0, 0.0, initialM_);//= beam + target;

  //final particles:
  unsigned int nParts = finalM_.size();
  Double_t masses[nParts];
  for(unsigned int i=0; i<nParts; i++){
    masses[i]=finalM_[i];
  }

  TGenPhaseSpace event;
  event.SetDecay(W, nParts, masses);

  TLorentzVector[nParts] *pFinal;
  TParticle fFinal[nParts];
  double weight, m23sq, m13sq, m12sq;
  //cout << "Start generation of y pi0 pi0 Dalitz" << endl;
  do{
      weight = event.Generate();


      for(unsigned int i=0; i<nParts; i++){
          pFinal[i] = event.GetDecay(i);
          fFinal[i] = TParticle(pids_[i],1,0,0,0,0,pFinal[i],W);
      }

      pPm23 = *pPim + *pPip;
      pPm13 = *pGamma + *pPim;

      m23sq=pPm23.M2(); m13sq=pPm13.M2();

      m12sq=M*M-m13sq-m23sq;
      if(m12sq<0){
        //cout << tmpm12_sq << "\t" << M*M << "\t" << m13_sq << "\t" << m23_sq << endl;
        //continue;
        m12sq=0.0001;
      }

      //call physics module
      vector<double> x;
      x.push_back(m23sq);
      x.push_back(m13sq);
      x.push_back(m12sq);
      double AMPpdf = model->intensity(x, minPar);

      double test = rando.Uniform(0,10);

      if(test<(weight*AMPpdf)){
          evt++;
          for(unsigned int i=0; i<nParts; i++)
            new((*fEvt)[i]) TParticle(fFinal[i]);

        fTree.Fill();
      }
  }while(evt<nEvents);

  fTree.Print();
  fTree.Write();
  output.Close();

  //cout << "Done ..." << endl << endl;

}
