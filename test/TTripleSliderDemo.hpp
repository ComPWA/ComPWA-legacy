// ComPWA header files
//#include "Core/Parameter.hpp"
//#include "Core/ParameterList.hpp"

// Root header files
#include <TApplication.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include "TGButton.h"
#include "TRootEmbeddedCanvas.h"
#include "TGLayout.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGTextEntry.h"
#include "TGSlider.h"
#include "TStyle.h"
#include "TGLabel.h"
#include "TGDoubleSlider.h"
#include <TQObject.h>
#include <TGButtonGroup.h>

#include "TComplex.h"
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TArrow.h"

//class TGWindow;
//class TGMainFrame;

class TTripleSliderDemo {//: public TGMainFrame {

  //RQ_OBJECT("TTripleSliderDemo");

  private:
     TGMainFrame         *panel;
     TRootEmbeddedCanvas *fCanvas;
     TGLayoutHints       *fLcan;
     TF1                 *fFcn[7];
     TGHorizontalFrame   *fHframe0, *fHframe1, *fHframe2;
     TGLayoutHints       *fBly, *fBfly1, *fBfly2, *fBfly3;
     TGHSlider           *fHSl[7];
     TGDoubleHSlider     *fSlRange;
     TGTextEntry         *fTEnt[7];
     TGTextBuffer        *fTBuf[7];
     TGLabel             *fLab[8];
     TGCheckButton       *fCheck1, *fCheck2, *fCheck3, *fCheck4, *fCheck5, *fCheck6;
     TGraph              *fArSum;//, *fAr1, fAr2;

     static int gmax;

     static double mass;
     static bool rel;
     static bool kmat;

     static Double_t breakup(Double_t M, Double_t m1, Double_t m2);
     static Double_t BarrierP(Double_t q, Double_t q0, Int_t L, Double_t d);
     static Double_t Barrier(Double_t q, Int_t L, Double_t d);
     static Double_t BWGamma(Double_t m, Double_t m0, Double_t Gamma0, Double_t m1, Double_t m2, Int_t L, Double_t d);
     static TComplex BreitWigner(Double_t m, Double_t m0, Double_t Gamma0, Double_t md1=0.135, Double_t md2=0.135, Int_t L=1, Double_t d=5.067);
     static TComplex BreitWignerKM(Double_t m, Double_t m0, Double_t Gamma0, Double_t m1, Double_t Gamma1, Double_t A0, Double_t A1, Int_t L=0, Double_t d=5.067);
     static TComplex Sum2BWC(double *x, double *p);
     static double Sum2BW(double *x, double *p);
     static double PhaseDiffBW(double *x, double *p);
     static TComplex BWC(double *x, double *p);
     static double BW(double *x, double *p);
     static double PhaseBW(double *x, double *p);
     double maxm(TF1 *f1, TF1 *f2, TF1 *f3);
     static void drawLabels(double A1, double A2, double ph);
     static void confFcn(TF1 *f, int col=1, int style = 1, TString tity="I (a.u.)", TString titx="m [GeV/c^2]");
     static double maxv(double x1, double x2, double x3);
     static int findNearestIndex(double *x, double *y, double m, int max);

  public:
     TTripleSliderDemo(double* par, double *slmin,double *slmax,double *step, TString *name, unsigned int nPar, bool inrel);
     virtual ~TTripleSliderDemo();

     void CloseWindow();
     void DoText();
     void DoSlider();
     void HandleButtons();
     TComplex ArgandPlot(TGraph *g, double xmin, double xmax, double *p, int nstp = 100);


     //ClassDef(TTripleSliderDemo, 0)
};

int TTripleSliderDemo::gmax = 500;

double TTripleSliderDemo::mass = 0.140;
bool TTripleSliderDemo::rel=false;
bool TTripleSliderDemo::kmat=false;

Double_t TTripleSliderDemo::breakup(Double_t M, Double_t m1, Double_t m2)
{
    Double_t sum12 = m1+m2;
    Double_t dif12 = m1-m2;

    Double_t q = 0.;

    if ( M > sum12 ) q = sqrt( (M*M - sum12*sum12) * (M*M - dif12*dif12) )/(2.*M);

    return q;
}

Double_t TTripleSliderDemo::BarrierP(Double_t q, Double_t q0, Int_t L, Double_t d)
{
    Double_t result = 1.0;

    Double_t z0 = q0*q0*d*d;
    Double_t z  = q*q*d*d;

    switch(L)
    {
    case 1: result = sqrt((1.+z0)/(1.+z));
            break;

    case 2: result = sqrt( ( (z0-3.)*(z0-3.)+9.*z0 ) / ( (z-3.)*(z-3.)+9.*z ) );
            break;

    default: break;
    }

    return result;
}

Double_t TTripleSliderDemo::Barrier(Double_t q, Int_t L, Double_t d)
{
    Double_t result = 1.0;

    Double_t z  = q*q*d*d;

    switch(L)
    {
    case 1: result = sqrt(2.*z/(1.+z));
            break;

    case 2: result = sqrt( ( 13.*z*z ) / ( (z-3.)*(z-3.)+9.*z ) );
            break;

    default: break;
    }

    return result;
}

Double_t TTripleSliderDemo::BWGamma(Double_t m, Double_t m0, Double_t Gamma0, Double_t m1, Double_t m2, Int_t L, Double_t d)
{
    Double_t q0 = breakup(m0,m1,m2);
    Double_t q  = breakup(m,m1,m2);

    if (0. == q) return 0.;

    Double_t phsp = q/q0;

    switch (L)
    {
    case 1: phsp = phsp*phsp*phsp;
            break;
    case 2: phsp = phsp*phsp*phsp*phsp*phsp;
            break;
    default: break;
    }

    Double_t bar = BarrierP(q,q0,L,d);

    Double_t G = Gamma0 * phsp * (m0/m) * bar*bar;

    return G;
}

TComplex TTripleSliderDemo::BreitWigner(Double_t m, Double_t m0, Double_t Gamma0, Double_t md1, Double_t md2, Int_t L, Double_t d) // d= 1/(0.197 GeV)
{
    md1=mass;
    md2=mass;

    TComplex I(0.,1.);
    TComplex BW, G;

    if (rel)
    {
        G=TComplex(BWGamma(m,m0,Gamma0,md1,md2,L,d),0);
        BW = m0*G/(m0*m0 - m*m - I*m0*G);
    }
    else
    {
        BW = Gamma0/2./(m0 - m - I*Gamma0/2.);
    }
    //if (damp) BW=breakup(m,m1,m2)*BW;

    return BW;
}

TComplex TTripleSliderDemo::BreitWignerKM(Double_t m, Double_t m0, Double_t Gamma0, Double_t m1, Double_t Gamma1, Double_t A0, Double_t A1, Int_t L, Double_t d) // d= 1/(0.197 GeV)
{

    TComplex I(0.,1.);

    double a0 = A0/sqrt(A0*A0+A1*A1);
    double a1 = A1/sqrt(A0*A0+A1*A1);

    TComplex G0(a0*BWGamma(m,m0,Gamma0,mass,mass,L,d),0);
    TComplex G1(a1*BWGamma(m,m1,Gamma1,mass,mass,L,d),0);

    double aux = (m0*m0-m*m)/(m1*m1-m*m);

    TComplex BW = m0*G0/(m0*m0 - m*m - I*m0*G0 - I*aux*m1*G1);

    return BW;
}

TComplex TTripleSliderDemo::Sum2BWC(double *x, double *p)
{
    double m = x[0];

    if (m<=0) return 0;

    double A0 = p[0];
    double m0 = p[1];
    double G0 = p[2];

    double A1 = p[3];
    double m1 = p[4];
    double G1 = p[5];

    double rp = p[6];

    double coh = p[7];

    TComplex rt(1,rp,true);

    TComplex bw0,bw1;

    if (kmat)
    {
        bw0 = BreitWignerKM(m, m0, G0, m1, G1, A0, A1);
        bw1 = BreitWignerKM(m, m1, G1, m0, G0, A0, A1);
    }
    else
    {
        bw0 = A0*BreitWigner(m, m0, G0);
        bw1 = A1*BreitWigner(m, m1, G1);
    }

    if (coh<0) return TComplex(sqrt(bw0.Rho2()+bw1.Rho2()),0);

    TComplex res = bw0+bw1*rt;

    return res;
}

double TTripleSliderDemo::Sum2BW(double *x, double *p)
{
    TComplex res = Sum2BWC(x,p);
    return res.Rho2();

    //double m = x[0];

    //if (m<=0) return 0;

    //double A0 = p[0];
    //double m0 = p[1];
    //double G0 = p[2];

    //double A1 = p[3];
    //double m1 = p[4];
    //double G1 = p[5];

    //double rp = p[6];

    //double coh = p[7];

    //TComplex rt(1,rp,true);

    //TComplex bw0 = A0*BreitWigner(m, m0, G0);
    //TComplex bw1 = A1*BreitWigner(m, m1, G1);

    //if (coh<0) return bw0.Rho2()+bw1.Rho2();

    //TComplex res = bw0+bw1*rt;

    //return res.Rho2();
}

double TTripleSliderDemo::PhaseDiffBW(double *x, double *p)
{
    double m = x[0];

    double A0 = p[0];
    double m0 = p[1];
    double G0 = p[2];

    double A1 = p[3];
    double m1 = p[4];
    double G1 = p[5];

    double rp = p[6];

    TComplex rt(1,rp,true);

    TComplex bw0,bw1;

    if (kmat)
    {
        bw0 = BreitWignerKM(m, m0, G0, m1, G1, A0, A1);
        bw1 = BreitWignerKM(m, m1, G1, m0, G0, A0, A1)*rt;
    }
    else
    {
        bw0 = A0*BreitWigner(m, m0, G0);
        bw1 = A1*BreitWigner(m, m1, G1)*rt;
    }

    double th0 = bw0.Theta();
    double th1 = bw1.Theta();

    double res = th0-th1;
    while (res>3.1415692) res -= 2*3.1415692;

    return fabs(res);
}

TComplex TTripleSliderDemo::BWC(double *x, double *p)
{
    double m = x[0];
    if (m<=0) return 0;

    double A0 = p[0];
    double m0 = p[1];
    double G0 = p[2];

    TComplex bw0 = A0*BreitWigner(m, m0, G0);

    return bw0;
}

double TTripleSliderDemo::BW(double *x, double *p)
{
    TComplex bw0 = BWC(x,p);
    return bw0.Rho2();
}



double TTripleSliderDemo::PhaseBW(double *x, double *p)
{
    double m = x[0];

    double A0 = p[0];
    double m0 = p[1];
    double G0 = p[2];

    TComplex bw0 = A0*BreitWigner(m, m0, G0);

    return bw0.Theta();
}

double TTripleSliderDemo::maxm(TF1 *f1, TF1 *f2, TF1 *f3)
{
    double x1 = f1->GetMaximum();
    double x2 = f2->GetMaximum();
    double x3 = f3->GetMaximum();

    if (x1>x2 && x1> x3) return x1;
    if (x2>x3) return x2;
    return x3;
}

void TTripleSliderDemo::drawLabels(double A1, double A2, double ph)
{
    TLatex lat;
    lat.SetNDC(kTRUE);

    double mt = gPad->GetTopMargin();
    double mb = gPad->GetBottomMargin();
    double ml = gPad->GetLeftMargin();
    double mr = gPad->GetRightMargin();

    lat.SetTextColor(1);
    lat.DrawLatex(ml+0.05*(1-mr-ml),0.97-mt-0.05, Form("A_{1} = %5.3f",A1));

    lat.SetTextColor(4);
    lat.DrawLatex(ml+0.05*(1-mr-ml),0.97-mt-0.11,  Form("A_{2} = %5.3f",A2));

    lat.SetTextColor(2);
    lat.DrawLatex(ml+0.05*(1-mr-ml),0.97-mt-0.17, Form("#Delta#phi = %5.1f#circ",ph));

}

void TTripleSliderDemo::confFcn(TF1 *f, int col, int style, TString tity, TString titx)
{
    f->SetNpx(1000);
    //f->GetHistogram()->SetXTitle(titx);
    //f->GetHistogram()->SetYTitle(tity);
    f->SetLineColor(col);
    f->SetLineStyle(style);
}

double TTripleSliderDemo::maxv(double x1, double x2, double x3)
{
    if (x1>x2 && x1>x3) return x1;
    if (x2>x1 && x2>x3) return x2;
    return x3;
}


int TTripleSliderDemo::findNearestIndex(double *x, double *y, double m, int max)
{
    double dx=1e8;
    int bidx=0;

    for (int i=0;i<max;++i) if (dx>fabs(x[i]-m)) {dx = fabs(x[i]-m); bidx=i;}
    return bidx;
}

TComplex TTripleSliderDemo::ArgandPlot(TGraph *g, double xmin, double xmax, double *p, int nstp)
{
    double dm = (xmax-xmin)/(nstp-1);

    int cnt = 0;
    TComplex res;

    for (int i=0;i<nstp;++i)
    {
        double m = xmin+i*dm;
        res = Sum2BWC(&m,p);
        g->SetPoint(cnt++, res.Re(), res.Im());
    }
    return res;
}

//______________________________________________________________________________
TTripleSliderDemo::TTripleSliderDemo(double *par,double *slmin,double *slmax,double *step, TString* name, unsigned int nPar, bool inrel)
 // : TGMainFrame(gClient->GetRoot(), 100, 100)
{
    panel = new TGMainFrame(gClient->GetRoot(),200,200);

    rel = inrel;
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    char buf[32];
    panel->SetCleanup(kDeepCleanup);
    // Create an embedded canvas and add to the main frame, centered in x and y
    // and with 30 pixel margins all around
    fCanvas = new TRootEmbeddedCanvas("Canvas", panel, 1200, 550);
    fLcan = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 10);
    panel->AddFrame(fCanvas, fLcan);
    //fCanvas->GetCanvas()->SetFillColor(33);
    //fCanvas->GetCanvas()->SetFrameFillColor(41);
    fCanvas->GetCanvas()->SetBorderMode(0);
    fCanvas->GetCanvas()->Divide(2,1);
    //fCanvas->GetCanvas()->GetPad(1)->SetLogy();
    fCanvas->GetCanvas()->GetPad(1)->SetGrid();
    fCanvas->GetCanvas()->GetPad(2)->SetGrid();

    //fCanvas->GetCanvas()->cd(1);

    int cnt=0;

    fHframe0 = new TGHorizontalFrame(panel, 0, 0, 0);
    fHframe1 = new TGHorizontalFrame(panel, 0, 0, 0);

    for(unsigned int ibut=0; ibut<nPar; ibut++){

      if(ibut<3){
        fLab[ibut]  = new TGLabel(fHframe0, name[ibut]);  fLab[ibut]->Resize(30,10);
        fHSl[ibut]  = new TGHSlider(fHframe0,200);        fHSl[ibut]->SetRange(slmin[ibut], slmax[ibut]);
        fTEnt[ibut] = new TGTextEntry(fHframe0, fTBuf[ibut] = new TGTextBuffer(8), cnt++);
        fTEnt[ibut]->Resize(50);
      }else{
        fLab[ibut]  = new TGLabel(fHframe1, name[ibut]);  fLab[ibut]->Resize(30,10);
        fHSl[ibut]  = new TGHSlider(fHframe1,200);        fHSl[ibut]->SetRange(slmin[ibut], slmax[ibut]);
        fTEnt[ibut] = new TGTextEntry(fHframe1, fTBuf[ibut] = new TGTextBuffer(8), cnt++);
        fTEnt[ibut]->Resize(50);
      }
    }


   /* fLab[0]  = new TGLabel(fHframe0, name[0]); fLab[0]->Resize(30,10);
    fHSl[0]  = new TGHSlider(fHframe0,200); fHSl[0]->SetRange(slmin[0], slmax[0]);
    fTEnt[0] = new TGTextEntry(fHframe0, fTBuf[0] = new TGTextBuffer(8), cnt++);
    fTEnt[0]->Resize(50);

    fLab[1] = new TGLabel(fHframe0, name[1]); fLab[1]->Resize(30,10);
    fHSl[1] = new TGHSlider(fHframe0,200); fHSl[1]->SetRange(slmin[1], slmax[1]);
    fTEnt[1] = new TGTextEntry(fHframe0, fTBuf[1] = new TGTextBuffer(8), cnt++);
    fTEnt[1]->Resize(50);

    fLab[2] = new TGLabel(fHframe0, name[2]); fLab[2]->Resize(30,10);
    fHSl[2] = new TGHSlider(fHframe0,200); fHSl[2]->SetRange(slmin[2], slmax[2]);
    fTEnt[2] = new TGTextEntry(fHframe0, fTBuf[2] = new TGTextBuffer(8), cnt++);
    fTEnt[2]->Resize(50);*/

    fLab[7]  = new TGLabel(fHframe0, "Range");
    fLab[7]->Resize(30,10);
    fSlRange = new TGDoubleHSlider(fHframe0, 100, kDoubleScaleBoth, cnt++);
    fSlRange->SetRange(0, 5);
    fSlRange->SetPosition(0, 5);
    fSlRange->Resize(215,100);
    fSlRange->Connect("PositionChanged()", "TTripleSliderDemo", panel, "DoSlider()");

    fCheck1 = new TGCheckButton(fHframe0, "Incoh", cnt++);
    fCheck1->SetState(kButtonUp);

    fCheck3 = new TGCheckButton(fHframe0, "Argand", cnt++);
    fCheck3->SetState(kButtonUp);

    fCheck6 = new TGCheckButton(fHframe0, "Vec", cnt++);
    fCheck6->SetState(kButtonUp);
    //fHframe0->Resize(200, 50);

    // ------------------------------------


    /*fLab[3] = new TGLabel(fHframe1, name[3]); fLab[3]->Resize(30,10);
    fHSl[3] = new TGHSlider(fHframe1,200); fHSl[3]->SetRange(slmin[3],slmax[3]);
    fTEnt[3] = new TGTextEntry(fHframe1, fTBuf[3] = new TGTextBuffer(8), cnt++);
    fTEnt[3]->Resize(50);

    fLab[4] = new TGLabel(fHframe1, name[4]); fLab[4]->Resize(30,10);
    fHSl[4] = new TGHSlider(fHframe1,200); fHSl[4]->SetRange(slmin[4],slmax[4]);
    fTEnt[4] = new TGTextEntry(fHframe1, fTBuf[4] = new TGTextBuffer(8), cnt++);
    fTEnt[4]->Resize(50);

    fLab[5] = new TGLabel(fHframe1, name[5]); fLab[5]->Resize(30,10);
    fHSl[5] = new TGHSlider(fHframe1,200); fHSl[5]->SetRange(slmin[5],slmax[5]);
    fTEnt[5] = new TGTextEntry(fHframe1, fTBuf[5] = new TGTextBuffer(8), cnt++);
    fTEnt[5]->Resize(50);

    fLab[6] = new TGLabel(fHframe1, name[6]); fLab[6]->Resize(30,10);
    fHSl[6] = new TGHSlider(fHframe1,200); fHSl[6]->SetRange(slmin[6],slmax[6]);
    fLab[6] = new TGLabel(fHframe1, "phase"); fLab[6]->Resize(30,10);
    fTEnt[6] = new TGTextEntry(fHframe1, fTBuf[6] = new TGTextBuffer(8), cnt++);
    fTEnt[6]->Resize(50);*/

    fCheck2 = new TGCheckButton(fHframe1, "Reson.", cnt++);
    fCheck2->SetState(kButtonUp);

    fCheck4 = new TGCheckButton(fHframe1, "rel", cnt++);
    fCheck4->SetState(rel?kButtonDown:kButtonUp);

    fCheck5 = new TGCheckButton(fHframe1, "KM", cnt++);
    fCheck5->SetState(kmat?kButtonDown:kButtonUp);

    fCheck1->Connect("Clicked()", "TTripleSliderDemo", panel, "DoSlider()");
    fCheck2->Connect("Clicked()", "TTripleSliderDemo", panel, "DoSlider()");
    fCheck3->Connect("Clicked()", "TTripleSliderDemo", panel, "DoSlider()");
    fCheck4->Connect("Clicked()", "TTripleSliderDemo", panel, "DoSlider()");
    fCheck5->Connect("Clicked()", "TTripleSliderDemo", panel, "DoSlider()");
    fCheck6->Connect("Clicked()", "TTripleSliderDemo", panel, "DoSlider()");

    // initial values
    fHSl[0]->SetPosition(par[0]);
    fHSl[1]->SetPosition(par[1]);
    fHSl[2]->SetPosition(par[2]);

    fHSl[3]->SetPosition(par[3]);
    fHSl[4]->SetPosition(par[4]);
    fHSl[5]->SetPosition(par[5]);
    fHSl[6]->SetPosition(par[6]);


    for (int i=0;i<7;++i)
    {
        fHSl[i]->Connect("PositionChanged(Int_t)", "TTripleSliderDemo", this, "DoSlider()");
        fHSl[i]->Resize(140,step[i]);
        //fHSl[i]->SetScale(step[i]);
        fTEnt[i]->Connect("ReturnPressed()", "TTripleSliderDemo", this, "DoText()");
    }

    //fBly = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 5, 5);
    fBly = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);

    fHframe0->AddFrame(fLab[0], fBly);
    fHframe0->AddFrame(fHSl[0], fBly);
    fHframe0->AddFrame(fTEnt[0], fBly);

    fHframe0->AddFrame(fLab[1], fBly);
    fHframe0->AddFrame(fHSl[1], fBly);
    fHframe0->AddFrame(fTEnt[1], fBly);

    fHframe0->AddFrame(fLab[2], fBly);
    fHframe0->AddFrame(fHSl[2], fBly);
    fHframe0->AddFrame(fTEnt[2], fBly);

    fHframe0->AddFrame(fLab[7], fBly);
    fHframe0->AddFrame(fSlRange, fBly);

    fHframe0->AddFrame(fCheck1, fBly);
    fHframe0->AddFrame(fCheck3, fBly);
    fHframe0->AddFrame(fCheck6, fBly);

    fHframe1->AddFrame(fLab[3], fBly);
    fHframe1->AddFrame(fHSl[3], fBly);
    fHframe1->AddFrame(fTEnt[3], fBly);

    fHframe1->AddFrame(fLab[4], fBly);
    fHframe1->AddFrame(fHSl[4], fBly);
    fHframe1->AddFrame(fTEnt[4], fBly);

    fHframe1->AddFrame(fLab[5], fBly);
    fHframe1->AddFrame(fHSl[5], fBly);
    fHframe1->AddFrame(fTEnt[5], fBly);

    fHframe1->AddFrame(fLab[6], fBly);
    fHframe1->AddFrame(fHSl[6], fBly);
    fHframe1->AddFrame(fTEnt[6], fBly);

    fHframe1->AddFrame(fCheck2, fBly);
    fHframe1->AddFrame(fCheck4, fBly);
    fHframe1->AddFrame(fCheck5, fBly);

    panel->AddFrame(fHframe0, fBly);
    panel->AddFrame(fHframe1, fBly);

    panel->SetWindowName("Interference Gui");
    panel->MapSubwindows();
    panel->Resize(panel->GetDefaultSize());
    panel->MapWindow();


    // ------------------------- The functions

    double mmin = 0;
    double mmax = 10;

    // BW sum coh
    fFcn[0] = new TF1("fSumBW",TTripleSliderDemo::Sum2BW,mmin,mmax,8);
    confFcn(fFcn[0],2);
    fFcn[0]->SetMinimum(0);

    // BW sum incoh
    fFcn[6] = new TF1("fSumBWi",TTripleSliderDemo::Sum2BW,mmin,mmax,8);
    confFcn(fFcn[6],6,2);
    fFcn[6]->SetMinimum(0);

    // phase sum
    fFcn[1] = new TF1("fSumPh",TTripleSliderDemo::PhaseDiffBW,mmin,mmax,7);
    confFcn(fFcn[1],2,1,"phase [rad]");
    fFcn[1]->SetMinimum(0);fFcn[1]->SetMaximum(3.2);

    // BW1
    fFcn[2] = new TF1("fBW1",TTripleSliderDemo::BW,mmin,mmax,3);
    confFcn(fFcn[2],1,9);
    fFcn[2]->SetMinimum(0);

    // BW2
    fFcn[3] = new TF1("fBW2",TTripleSliderDemo::BW,mmin,mmax,3);
    confFcn(fFcn[3],4,7);
    fFcn[3]->SetMinimum(0);

    // ph1
    fFcn[4] = new TF1("fPh1",TTripleSliderDemo::PhaseBW,mmin,mmax,3);
    confFcn(fFcn[4],1,9,"phase [rad]");
    fFcn[4]->SetMinimum(0);fFcn[4]->SetMaximum(3.2);

    // ph2
    fFcn[5] = new TF1("fPh2",TTripleSliderDemo::PhaseBW,mmin,mmax,3);
    confFcn(fFcn[5],4,7,"phase [rad]");
    fFcn[5]->SetMinimum(0);fFcn[5]->SetMaximum(3.2);

    // ------------------------- The graphs

    fArSum=new TGraph(gmax);
    fArSum->SetLineWidth(2);
    fArSum->SetLineColor(2);
    fArSum->SetMarkerColor(2);
    fArSum->SetMarkerSize(0.5);
    fArSum->SetMarkerStyle(20);

    //fAr1=new TGraph(500);
    //fAr1->SetLineWidth(2);
    //fAr1->SetLineColor(1);
    //fAr1->SetLineStyle(2);

    //fAr2=new TGraph(500);
    //fAr2->SetLineWidth(2);
    //fAr2->SetLineColor(4);
    //fAr2->SetLineStyle(7);



    //double A1 = fHSl[0]->GetPosition()/1000.;
    //double m1 = fHSl[1]->GetPosition()/1000.;
    //double G1 = fHSl[2]->GetPosition()/1000.;

    //double A2 = fHSl[3]->GetPosition()/1000.;
    //double m2 = fHSl[4]->GetPosition()/1000.;
    //double G2 = fHSl[5]->GetPosition()/1000.;

    //double ph = fHSl[6]->GetPosition();


    //fFcn[0]->SetParameters(A1, m1, G1, A2, m2, G2, ph/57.3);
    //fFcn[0]->SetRange(mmin,mmax);

    //fFcn[1]->SetParameters(A1, m1, G1, A2, m2, G2, ph/57.3);
    //fFcn[1]->SetRange(mmin,mmax);

    //fFcn[2]->SetParameters(A1, m1, G1); fFcn[2]->SetRange(mmin,mmax);
    //fFcn[3]->SetParameters(A2, m2, G2); fFcn[3]->SetRange(mmin,mmax);

    //fFcn[4]->SetParameters(A1, m1, G1); fFcn[4]->SetRange(mmin,mmax);
    //fFcn[5]->SetParameters(A2, m2, G2); fFcn[5]->SetRange(mmin,mmax);

    TTripleSliderDemo::DoSlider();

//  fHframe0->Resize(200, 50);


/*
    fCheck1 = new TGCheckButton(fHframe0, "&Constrained", HCId1);
    fCheck2 = new TGCheckButton(fHframe0, "&Relative", HCId2);
    fCheck1->SetState(kButtonUp);
    fCheck2->SetState(kButtonUp);
    fCheck1->SetToolTipText("Pointer position constrained to slider sides");
    fCheck2->SetToolTipText("Pointer position relative to slider position");



    fHframe0->Resize(200, 50);

    fHframe1 = new TGHorizontalFrame(this, 0, 0, 0);

    fHslider1 = new TGTripleHSlider(fHframe1, 190, kDoubleScaleBoth, HSId1,
                                   kHorizontalFrame,
                                   GetDefaultFrameBackground(),
                                   kFALSE, kFALSE, kFALSE, kFALSE);
    fHslider1->Connect("PointerPositionChanged()", "TTripleSliderDemo",
                      this, "DoSlider()");
    fHslider1->Connect("PositionChanged()", "TTripleSliderDemo",
                      this, "DoSlider()");
    fHslider1->SetRange(0.05,5.0);

    fHframe1->Resize(200, 25);

    fHframe2 = new TGHorizontalFrame(this, 0, 0, 0);

    fTeh1 = new TGTextEntry(fHframe2, fTbh1 = new TGTextBuffer(5), HId1);
    fTeh2 = new TGTextEntry(fHframe2, fTbh2 = new TGTextBuffer(5), HId2);
    fTeh3 = new TGTextEntry(fHframe2, fTbh3 = new TGTextBuffer(5), HId3);

    fTeh1->SetToolTipText("Minimum (left) Value of Slider");
    fTeh2->SetToolTipText("Pointer Position Value");
    fTeh3->SetToolTipText("Maximum (right) Value of Slider");

    fTbh1->AddText(0, "0.0");
    fTbh2->AddText(0, "0.0");
    fTbh3->AddText(0, "0.0");

    fTeh1->Connect("TextChanged(char*)", "TTripleSliderDemo", this,
                  "DoText(char*)");
    fTeh2->Connect("TextChanged(char*)", "TTripleSliderDemo", this,
                  "DoText(char*)");
    fTeh3->Connect("TextChanged(char*)", "TTripleSliderDemo", this,
                  "DoText(char*)");

    fCheck1->Connect("Clicked()", "TTripleSliderDemo", this,
                    "HandleButtons()");
    fCheck2->Connect("Clicked()", "TTripleSliderDemo", this,
                    "HandleButtons()");

    fHframe2->Resize(100, 25);

    //--- layout for buttons: top align, equally expand horizontally
    fBly = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 5, 5);

    //--- layout for the frame: place at bottom, right aligned
    fBfly1 = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 5, 5, 5, 5);
    fBfly2 = new TGLayoutHints(kLHintsTop | kLHintsLeft,    5, 5, 5, 5);
    fBfly3 = new TGLayoutHints(kLHintsTop | kLHintsRight,   5, 5, 5, 5);

    fHframe0->AddFrame(fCheck1, fBfly2);
    fHframe0->AddFrame(fCheck2, fBfly2);
    fHframe1->AddFrame(fHslider1, fBly);
    fHframe2->AddFrame(fTeh1, fBfly2);
    fHframe2->AddFrame(fTeh2, fBfly1);
    fHframe2->AddFrame(fTeh3, fBfly3);

    AddFrame(fHframe0, fBly);
    AddFrame(fHframe1, fBly);
    AddFrame(fHframe2, fBly);

    // Set main frame name, map sub windows (buttons), initialize layout
    // algorithm via Resize() and map main frame
    SetWindowName("Triple Slider Demo");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();

    fFitFcn = new TF1("fFitFcn", "TMath::LogNormal(x, [0], [1], [2])", 0, 5);
    fFitFcn->SetRange(0.0, 2.5);
    fFitFcn->SetParameters(1.0, 0, 1);
    fFitFcn->SetMinimum(1.0e-3);
    fFitFcn->SetMaximum(10.0);
    fFitFcn->SetLineColor(kRed);
    fFitFcn->SetLineWidth(1);
    fFitFcn->Draw();

    fHslider1->SetPosition(0.05,2.5);
    fHslider1->SetPointerPosition(1.0);

    sprintf(buf, "%.3f", fHslider1->GetMinPosition());
    fTbh1->Clear();
    fTbh1->AddText(0, buf);
    sprintf(buf, "%.3f", fHslider1->GetPointerPosition());
    fTbh2->Clear();
    fTbh2->AddText(0, buf);
    sprintf(buf, "%.3f", fHslider1->GetMaxPosition());
    fTbh3->Clear();
    fTbh3->AddText(0, buf);
    */
}

void TTripleSliderDemo::DoSlider()
{
   // Handle slider widgets.

    rel  = fCheck4->GetState();
    kmat = fCheck5->GetState();

    double A1 = fHSl[0]->GetPosition()/1000.;
    double m1 = fHSl[1]->GetPosition()/1000.;
    double G1 = fHSl[2]->GetPosition()/1000.;

    double A2 = fHSl[3]->GetPosition()/1000.;
    double m2 = fHSl[4]->GetPosition()/1000.;
    double G2 = fHSl[5]->GetPosition()/1000.;

    double ph = fHSl[6]->GetPosition();

    float mmin, mmax;

    fSlRange->GetPosition(mmin,mmax);

    //cout <<mmin<<" "<<mmax<<endl;

    for (int i=0;i<7;++i)
    {
        TString txt=Form("%d",fHSl[i]->GetPosition());
        fTBuf[i]->Clear();
        fTBuf[i]->AddText(0,txt.Data());
        gClient->NeedRedraw(fTEnt[i]);
    }


    if (mmin<0.01 && mmax>4.99)
    {
        mmin = m1-4*(G1+G2);
        mmax = m2+4*(G1+G2);
        if (m1>m2)
        {
            mmin = m2-4*(G1+G2);
            mmax = m1+4*(G1+G2);
        }
    }

    if (mmin<2*mass) mmin=2*mass;

    fFcn[0]->SetParameters(A1, m1, G1, A2, m2, G2, ph/57.3,1);
    fFcn[0]->SetRange(mmin,mmax);

    fFcn[6]->SetParameters(A1, m1, G1, A2, m2, G2, ph/57.3,-1);
    fFcn[6]->SetRange(mmin,mmax);

    fFcn[1]->SetParameters(A1, m1, G1, A2, m2, G2, ph/57.3);
    fFcn[1]->SetRange(mmin,mmax);

    fFcn[2]->SetParameters(A1, m1, G1); fFcn[2]->SetRange(mmin,mmax);
    fFcn[3]->SetParameters(A2, m2, G2); fFcn[3]->SetRange(mmin,mmax);

    fFcn[4]->SetParameters(A1, m1, G1); fFcn[4]->SetRange(mmin,mmax);
    fFcn[5]->SetParameters(A2, m2, G2); fFcn[5]->SetRange(mmin,mmax);

    fCanvas->GetCanvas()->cd(1);

    //double fmax = maxv(fFcn[0]->GetMaximum(), fFcn[2]->GetMaximum(), fFcn[3]->GetMaximum()) ;

    double fmax = maxm(fFcn[0], fFcn[6], fFcn[3])*1.05;

    fFcn[6]->SetMaximum(fmax);
    fFcn[0]->SetMaximum(fmax);

    fFcn[6]->GetHistogram()->SetXTitle("m [GeV/c^{2}]");
    fFcn[6]->GetHistogram()->SetYTitle("intensity [A.U.]");
    fFcn[0]->GetHistogram()->SetXTitle("m [GeV/c^{2}]");
    fFcn[0]->GetHistogram()->SetYTitle("intensity [A.U.]");

    fFcn[1]->GetHistogram()->SetXTitle("m [GeV/c^{2}]");
    fFcn[1]->GetHistogram()->SetYTitle("phase [rad]");


    TString opt = "";
    if (fCheck1->GetState())
    {
        fFcn[6]->Draw();
        opt="same";
    }

    fFcn[0]->Draw(opt);
    opt="same";

    if (fCheck2->GetState())
    {
        fFcn[2]->Draw(opt);
        fFcn[3]->Draw(opt);
    }

    TLegend leg(0.7, 0.8, 1-gPad->GetRightMargin(), 1- gPad->GetTopMargin());
    leg.AddEntry(fFcn[0],"|A1+A2|^2");
    if (fCheck1->GetState())  leg.AddEntry(fFcn[6],"I1 + I2");
    if (fCheck2->GetState()) {leg.AddEntry(fFcn[2],"I1");leg.AddEntry(fFcn[3],"I2");}

    if (fCheck2->GetState()||fCheck1->GetState()) leg.Draw();
    fCanvas->GetCanvas()->Update();

    fCanvas->GetCanvas()->cd(2);

    if (fCheck3->GetState())
    {
        TTripleSliderDemo::ArgandPlot(fArSum, mmin, mmax, fFcn[0]->GetParameters(), gmax);
        fArSum->GetHistogram()->SetTitle("Argand Diagram");
        fArSum->GetHistogram()->SetXTitle("Re(A1+A2)");
        fArSum->GetHistogram()->SetYTitle("Im(A1+A2)");

        fArSum->Draw("APC");

        TMarker ma;
        ma.SetMarkerSize(1.5);
        ma.SetMarkerStyle(21);
        ma.SetMarkerColor(2);

        double *gx=fArSum->GetX();
        double *gy=fArSum->GetY();

        ma.DrawMarker(gx[0],gy[0]);

        int idx1 = (double)gmax*(m1-mmin)/(mmax-mmin);//findNearestIndex(x, y, m1, gmax);
        int idx2 = (double)gmax*(m2-mmin)/(mmax-mmin);//findNearestIndex(x, y, m2, gmax);

        ma.SetMarkerStyle(20);
        ma.SetMarkerColor(1);
        ma.DrawMarker(gx[idx1],gy[idx1]);

        ma.SetMarkerColor(4);
        ma.DrawMarker(gx[idx2],gy[idx2]);

        if (fCheck6->GetState())
        {
            // some vector painting

            TArrow ar1, ar2, ar3;

            ar1.SetLineColor(1);
            ar2.SetLineColor(4);
            ar3.SetLineColor(2);
            ar1.SetFillColor(1);
            ar2.SetFillColor(4);
            ar3.SetFillColor(2);
            ar3.SetLineWidth(2);

            const int nstp = 10;
            double dm = (mmax-mmin)/(nstp-1);
            double asz=.007;
            double offy = fmax*0.5;
            TComplex rt(1,ph/57.3,true);
           // cout.precision(3);

            double xsf = 1.;///fmax;
            double ysf = 1.;///(mmax-mmin);

            double mar[nstp+2];

            for (int i=0;i<nstp;++i) mar[i] = mmin+i*dm;
            mar[nstp]   = m1;
            mar[nstp+1] = m2;

            for (int i=0;i<nstp+2;++i)
            {
                double m = mar[i];
                TComplex bw0 = BWC(&m,fFcn[2]->GetParameters());
                TComplex bw1 = BWC(&m,fFcn[3]->GetParameters())*rt;
                TComplex tot = Sum2BWC(&m,fFcn[0]->GetParameters());

                double x0 = bw0.Re()*xsf;
                double y0 = bw0.Im()*ysf;

                double x1 = bw1.Re()*xsf;
                double y1 = bw1.Im()*ysf;

                double xt = tot.Re()*xsf;
                double yt = tot.Im()*ysf;

                //if (i==5)
                //{
                    ////cout <<"("<<i<<": "<<bw1.Re()<<","<<bw1.Im()<<" / "<<bw1.Rho()<<","<<bw1.Theta()<<") * ("<<rt.Re()<<","<<rt.Im()<<" / "<<rt.Rho()<<","<<rt.Theta()<<")";
                    ////bw1 = bw1 * rt;
                    ////cout <<" --> ("<<x1<<","<<y1<<" / "<<bw1.Rho()<<","<<bw1.Theta()<<")"<<endl;
                //ar1.DrawArrow(m,    offy,    m+x1,    offy+y1,     asz);
                //}

                //ar1.DrawArrow(m,    offy,    m+x0,    offy+y0,     asz);

                //ar2.DrawArrow(m+x0,    offy+y0, m+x1+x0, offy+y1+y0,  asz);

                //ar3.DrawArrow(m,    offy,    m+xt,    offy+yt,     asz);

                ar1.DrawArrow(0,   0,    x0,    y0,     asz);
                ar2.DrawArrow(x0, y0, x1+x0, y1+y0,  asz);
                //ar3.DrawArrow(0,  0 ,    xt,    yt,     asz);

                //ar3.DrawArrow(m,    offy, m+x0+x1  ,   offy+y0+y1 , asz);
            }
        }

        //ArgantPlot(fAr1, mmin, mmax, fFcn[0]->GetParameters(), 100);
        //fAr1->Draw("C same");

        //ArgantPlot(fAr, mmin, mmax, fFcn[0]->GetParameters(), 100);
        //fAr2->Draw("C same");
    }
    else
    {
        fFcn[1]->Draw();
        fFcn[4]->Draw("same");
        fFcn[5]->Draw("same");

    }


    fCanvas->GetCanvas()->cd();

    fCanvas->GetCanvas()->Update();



   //char buf[32];

   //sprintf(buf, "%.3f", fHslider1->GetMinPosition());
   //fTbh1->Clear();
   //fTbh1->AddText(0, buf);
   //fTeh1->SetCursorPosition(fTeh1->GetCursorPosition());
   //fTeh1->Deselect();
   //gClient->NeedRedraw(fTeh1);

   //sprintf(buf, "%.3f", fHslider1->GetPointerPosition());
   //fTbh2->Clear();
   //fTbh2->AddText(0, buf);
   //fTeh2->SetCursorPosition(fTeh2->GetCursorPosition());
   //fTeh2->Deselect();
   //gClient->NeedRedraw(fTeh2);

   //sprintf(buf, "%.3f", fHslider1->GetMaxPosition());
   //fTbh3->Clear();
   //fTbh3->AddText(0, buf);
   //fTeh3->SetCursorPosition(fTeh3->GetCursorPosition());
   //fTeh3->Deselect();
   //gClient->NeedRedraw(fTeh3);

   //fFitFcn->SetParameters(fHslider1->GetPointerPosition(), 0, 1);
   //fFitFcn->SetRange(fHslider1->GetMinPosition()-0.05,
                     //fHslider1->GetMaxPosition());
   //fFitFcn->Draw();
   //fCanvas->GetCanvas()->Modified();
   //fCanvas->GetCanvas()->Update();
}

//______________________________________________________________________________
void TTripleSliderDemo::HandleButtons()
{
   // Handle different buttons.

   //TGButton *btn = (TGButton *) gTQSender;
   //Int_t id = btn->WidgetId();

   //switch (id) {
      //case HCId1:
         //fHslider1->SetConstrained(fCheck1->GetState());
         //break;
      //case HCId2:
         //fHslider1->SetRelative(fCheck2->GetState());
         //break;
      //default:
         //break;
   //}
}

//______________________________________________________________________________
 TTripleSliderDemo::~TTripleSliderDemo()
 {
    // Clean up

   panel->Cleanup();
 }

 //______________________________________________________________________________
 void TTripleSliderDemo::CloseWindow()
 {
    // Called when window is closed via the window manager.

    delete panel;
 }

 //______________________________________________________________________________
 void TTripleSliderDemo::DoText()
 {
    // Handle text entry widgets.

     for (int i=0;i<7;++i){
         fHSl[i]->SetPosition(Int_t(atof(fTBuf[i]->GetString())));
         //cout <<Int_t(atof(fTBuf[i]->GetString()))<<endl;
     }
     TTripleSliderDemo::DoSlider();


    //TGTextEntry *te = (TGTextEntry *) gTQSender;
    //Int_t id = te->WidgetId();

    //switch (id) {
       //case HId1:
          //fHslider1->SetPosition(atof(fTbh1->GetString()),
                                 //fHslider1->GetMaxPosition());
          //break;
       //case HId2:
          //fHslider1->SetPointerPosition(atof(fTbh2->GetString()));
          //break;
       //case HId3:
          //fHslider1->SetPosition(fHslider1->GetMinPosition(),
                                 //atof(fTbh1->GetString()));
          //break;
       //default:
          //break;
    //}
    //fFitFcn->SetParameters(fHslider1->GetPointerPosition(), 0, 1);
    //fFitFcn->SetRange(fHslider1->GetMinPosition()-0.05,
                      //fHslider1->GetMaxPosition());

    //fCanvas->GetCanvas()->GetPad(1)->cd(1);//->SetLogy(logscale);
    //fFitFcn->Draw();

    //fCanvas->GetCanvas()->Modified();
    //fCanvas->GetCanvas()->Update();
 }
