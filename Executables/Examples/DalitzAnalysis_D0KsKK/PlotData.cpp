
#include <stdio.h>
#include <numeric>

#include "Core/DataPoint.hpp"
#include "Core/Logging.hpp"
#include "Core/ProgressBar.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"
#include "PlotData.hpp"

#include "Math/ProbFuncMathCore.h"
#include "TLegend.h"

#include "histTools.hpp"

using namespace ComPWA;

plotData::plotData(std::string name, int bins) : _name(name),
_bins(bins), _isFilled(0), _globalScale(1.0),
dataDiagrams("data","Data", bins),
phspDiagrams("phsp","Phase-space",bins),
fitDiagrams("fit","Model",bins),
fitHitMissDiagrams("fitHitMiss","HitMiss",bins),
h_weights("h_weights","h_weights",bins,0,1.01)
{
	gStyle->SetOptStat(10);//entries only
	//	gStyle->SetOptStat(1000001); //name and integral
	gStyle->SetOptTitle(0);

	//Full intensity blue
	fitDiagrams.setColor(kBlue-4);
	//Phase-space green
	phspDiagrams.setColor(kGreen);

	fitDiagrams.SetStats(0);
	phspDiagrams.SetStats(0);


	ComPWA::Physics::DPKinematics::DalitzKinematics* kin =
			dynamic_cast<ComPWA::Physics::DPKinematics::DalitzKinematics*>(Kinematics::instance());

	//=== generate contour
	double xpoints[4001],ypoints[4001];
	kin->phspContour(0,1,2000,xpoints,ypoints);
	m23m13_contour = TGraph(4001,xpoints,ypoints);
	m23m13_contour.SetMarkerStyle(1);
	m23m13_contour.SetLineColor(kRed);
	m23m13_contour.SetMarkerColor(kRed);
	m23m13_contour.SetMarkerSize(0.0);
	m23m13_contour.SetTitle("phspContour");
	m23m13_contour.SetFillColor(kWhite);
	kin->phspContour(0,2,2000,xpoints,ypoints);
	m23m12_contour = TGraph(4001,xpoints,ypoints);
	m23m12_contour.SetMarkerStyle(1);
	m23m12_contour.SetLineColor(kRed);
	m23m12_contour.SetMarkerColor(kRed);
	m23m12_contour.SetMarkerSize(0.0);
	m23m12_contour.SetTitle("phspContour");
	m23m12_contour.SetFillColor(kWhite);
	kin->phspContour(2,1,2000,xpoints,ypoints);
	m12m13_contour = TGraph(4001,xpoints,ypoints);
	m12m13_contour.SetMarkerStyle(1);
	m12m13_contour.SetLineColor(kRed);
	m12m13_contour.SetMarkerColor(kRed);
	m12m13_contour.SetMarkerSize(0.0);
	m12m13_contour.SetTitle("phspContour");
	m12m13_contour.SetFillColor(kWhite);

}

plotData::~plotData()
{

}

TH2Poly* plotData::getAdBinHist(int bins)
{
	BOOST_LOG_TRIVIAL(info)<<"Creating adaptively binned histogram";

	if( !s_data )
		throw std::runtime_error("plotData::getAdBinHist() | No data sample set!");

	//Which sample do we use for calculation of binning? Use Hit&Miss sample of available
    std::shared_ptr<ComPWA::DataReader::Data> adaptiveBinningSample = s_data;
	//	if(s_hitMiss) adaptiveBinningSample = s_hitMiss;

	Double_t* tmpArray = new Double_t[2*(adaptiveBinningSample->getNEvents())]; //create array on the heap
	for(unsigned int i = 0; i < adaptiveBinningSample->getNEvents(); i++){//loop over data
		dataPoint point(adaptiveBinningSample->getEvent(i));
		double mSq = point.getVal(0);
		double phi = point.getVal(8);
		tmpArray[i]=mSq; //mKKsq
		tmpArray[adaptiveBinningSample->getNEvents()+i]=phi; //mKSK+sq
	}
	//Create adaptive histogram for data
	TH2Poly* hist = adaptiveBinning(adaptiveBinningSample->getNEvents(),2,tmpArray,bins);
	if( !hist )
		throw std::runtime_error("plotData::getAdBinHist() | Can not create"
				"adaptively binned histogram!");

	hist->SetName("dataAdaptiveBinned");
	hist->GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}]");
	hist->GetXaxis()->SetNdivisions(508);
	hist->GetXaxis()->SetTitleSize(0.06);
	hist->GetXaxis()->SetLabelSize(0.05);
	hist->GetYaxis()->SetTitle("cos#Theta_{KK}");
	hist->GetYaxis()->SetRangeUser(-1,1);
	hist->GetYaxis()->SetTitleSize(0.06);
	hist->GetYaxis()->SetLabelSize(0.05);
	hist->Reset("");

	return hist;
}

void plotData::setFitAmp(std::vector<std::shared_ptr<Amplitude> > ampVec,
		std::vector<double> fraction)
{
	_ampVec=ampVec;
	_fraction=fraction;

	double sumFraction = std::accumulate(_fraction.begin(), _fraction.end(), 0.0);

	if(sumFraction > 1.0)
		throw std::runtime_error("plotData::setFitAmp() | Fractions sum larger 1.0!");

	if(_fraction.size() == _ampVec.size()-1)
		_fraction.push_back(1-sumFraction);

	//check size
	if( _fraction.size() != _ampVec.size() )
		throw std::runtime_error("plotData::setFitAmp() | List of fractions "
				"does not match with list of amplitudes!");


	auto it = _ampVec.begin();
	for( ; it!=_ampVec.end(); ++it){
		ampHistos.push_back(
				dalitzHisto( (*it)->GetName(), (*it)->GetName(), _bins )
		);
		ampHistos.back().setColor(kBlack);
		ampHistos.back().SetStats(0);
	}
	resonanceItr resit = _ampVec.at(0)->GetResonanceItrFirst();
	for( ; resit!=_ampVec.at(0)->GetResonanceItrLast(); ++resit){
		if( (*resit)->GetName().find("_CP") != std::string::npos ) continue;
		signalComponents.push_back(
				dalitzHisto( (*resit)->GetName(), (*resit)->GetName(), _bins )
		);
		signalComponents.back().setColor(kBlack);
		signalComponents.back().SetStats(0);
	}
	_isFilled = 0;
}

void plotData::Fill()
{
	//TODO: reset diagrams here
	try{
		dataAdaptiveBinningHist = getAdBinHist(_bins);
		fitAdaptiveBinningHist =
				(TH2Poly*) dataAdaptiveBinningHist->Clone("fitAdaptiveBinned");
	} catch (std::exception& ex){
		BOOST_LOG_TRIVIAL(error) <<"plotData::Fill() | "<<ex.what();
		throw;
	}

	//===== Fill data histograms
	if( s_data ){
		for(unsigned int i = 0; i < s_data->getNEvents(); i++){//loop over data
			Event event(s_data->getEvent(i));

			double eff = 1.0;
			if( _correctForEfficiency ) eff = event.getEfficiency();
			if( eff==0.0 ){
				BOOST_LOG_TRIVIAL(error) <<"plotData::Fill() | Loop over "
						"data sample: An event with zero efficiency was found! "
						"This should not happen! We skip it!";
				continue;
			}
			double evWeight = event.getWeight();

			dataDiagrams.Fill(event,evWeight*1/eff);
			h_weights.Fill(evWeight*1/eff);

			dataPoint point(event);
			dataAdaptiveBinningHist->Fill(
					point.getVal(0),
					point.getVal(8),
					evWeight*1/eff
			);
		}
		_globalScale = dataDiagrams.GetIntegral();
	}

	//===== Plot amplitude
	BOOST_LOG_TRIVIAL(info) << "PlotData::plot | Evaluating amplitude...";

	if( _ampVec.size() && s_phsp ){
		std::vector<double> ampHistos_integral(_ampVec.size(),0);

		progressBar bar(s_phsp->getNEvents());
		for(unsigned int i = 0; i < s_phsp->getNEvents(); i++){//loop over phsp MC
			bar.nextEvent();
			Event event(s_phsp->getEvent(i));
			double eff = 1.0;
			if( _correctForEfficiency ) eff = event.getEfficiency();
			if( eff==0.0 ){
				BOOST_LOG_TRIVIAL(error) <<"plotData::Fill() | Loop over "
						"phsp sample: An event with zero efficiency was found! "
						"This should not happen! We skip it!";
				continue;
			}
			double evWeight = event.getWeight();

			phspDiagrams.Fill(event,evWeight*1/eff);//scale phsp to data size

			dataPoint point;
			try{
				point = dataPoint(event);
			} catch (BeyondPhsp& ex){ //event outside phase, remove
				continue;
			}

			for(int t=0; t< signalComponents.size(); ++t){
				resonanceItr res = _ampVec.at(0)->GetResonanceItrFirst();
				for(int j=0; j<t; ++j) res++;

				//skip CP partner of resonance
				if( (*res)->GetName().find("_CP") != std::string::npos )
					continue;

				std::complex<double> val = (*res)->Evaluate(point);

				try{ // trying to find a CP partner and add it
					std::shared_ptr<Resonance> cpRes;
					std::shared_ptr<ComPWA::Physics::AmplitudeSum::AmpSumIntensity> tmpAmp =
							std::dynamic_pointer_cast<ComPWA::Physics::AmplitudeSum::AmpSumIntensity>(
									_ampVec.at(0)
							);
					cpRes = tmpAmp->GetResonance(
							(*res)->GetName()+"_CP");
					val += cpRes->Evaluate(point);
				} catch(std::exception& ex){ }
				signalComponents.at(t).Fill( event, std::norm(val)/evWeight*1/eff );
			}

			std::complex<double> tmp_intens2(0,0);
			auto it = _ampVec.at(0)->GetResonanceItrFirst();
			for( ; it != _ampVec.at(0)->GetResonanceItrLast(); ++it){
				if( (*it)->GetName().find("_CP") != std::string::npos ) continue;
				tmp_intens2 += (*it)->Evaluate(point);
			}
			ampHistos.at(0).Fill( event, std::norm(tmp_intens2)/evWeight*1/eff );

			double intens = 0;
			for(int t=0; t<_ampVec.size(); ++t){
				ParameterList tmp_list = _ampVec.at(t)->intensity(point);
				double tmp_intens = tmp_list.GetDoubleParameter(0)->GetValue();
				tmp_intens = tmp_intens	*_fraction.at(t);
				intens += tmp_intens;
				if( t==1 ) //fill background
					ampHistos.at(t).Fill( event, tmp_intens/evWeight*1/eff );
			}
			fitDiagrams.Fill( event, intens/evWeight*1/eff );

			//Fill adaptive histogram
			fitAdaptiveBinningHist->Fill(
					point.getVal(0), point.getVal(8), intens/evWeight*1/eff
			);
		}

		_globalScale = 1934;
		//Scale histograms to match data sample
		fitDiagrams.Scale( _globalScale /fitDiagrams.GetIntegral() );
		phspDiagrams.Scale( _globalScale /phspDiagrams.GetIntegral() );

		for(int t=0; t<_ampVec.size(); ++t){
			double scale = _globalScale
					/ ampHistos.at(t).GetIntegral()
					* _fraction.at(t);
			ampHistos.at(t).Scale( scale );
		}
		for(int t=0; t<signalComponents.size(); ++t){
			double scale = _globalScale
					/ ampHistos.at(0).GetIntegral()
					* _fraction.at(0);
			signalComponents.at(t).Scale( scale );
		}
	}

	//===== Plot hit&miss data
	if(s_hitMiss){
		for(unsigned int i = 0; i < s_hitMiss->getNEvents(); i++){//loop over data
			Event event(s_hitMiss->getEvent(i));
			double eff = 1.0;
			if( _correctForEfficiency ) eff = event.getEfficiency();
			if( eff==0.0 ){
				BOOST_LOG_TRIVIAL(error) <<"plotData::Fill() | Loop over "
						"Hit&Miss sample: An event with zero efficiency was "
						"found! This should not happen! We skip it!";
				continue;
			}
			double evWeight = event.getWeight();

			fitHitMissDiagrams.Fill(event,evWeight*1/eff);
		}
	}

	_isFilled = 1;
}

void plotData::Plot()
{
	if( !_isFilled ) Fill();

	ComPWA::Physics::DPKinematics::DalitzKinematics* kin =
			dynamic_cast<ComPWA::Physics::DPKinematics::DalitzKinematics*>(Kinematics::instance());

	//=== create 2dim residual histogrom from data and fit model. Use adaptive binning here
	TH2Poly* adaptiveResiduals =
			(TH2Poly*) dalitzHisto::getTH2PolyPull(
					dataAdaptiveBinningHist,
					fitAdaptiveBinningHist
			);
	adaptiveResiduals->SetStats(0);

	//=== calculate goodness-of-fit variable
	double chi2;
	int NDF, igood;
	//	dalitzHisto::getTH2PolyChi2(dataAdaptiveBinningHist,fitAdaptiveBinningHist,chi2,NDF,igood);
	//	double pValue = ROOT::Math::chisquared_cdf_c(chi2,NDF);//lower tail cdf
	//	BOOST_LOG_TRIVIAL(info) << "Goodness-of-fit (2D): chi^2/ndf = "
	//			<<chi2<<" / "<<NDF<<" = "<<chi2/NDF<<" pValue="<<pValue;

	Chi2TestX(dataAdaptiveBinningHist, fitAdaptiveBinningHist,chi2,NDF,igood,"UW,UF");
	double pValue = ROOT::Math::chisquared_cdf_c(chi2,NDF);//lower tail cdf
	BOOST_LOG_TRIVIAL(info) << "Goodness-of-fit (2D): chi^2/ndf = "
			<<chi2<<" / "<<NDF<<" = "<<chi2/NDF<<" pValue="<<pValue;

	char chi2Char[60];
	sprintf(chi2Char,"#chi^{2}/ndf = %.2f/%d",chi2,NDF);


	//----- plotting dalitz distributions -----
	TCanvas* c1 = new TCanvas("dalitz","dalitz",50,50,800,800);
	c1->Divide(2,2);
	c1->cd(1);
	dataDiagrams.getHistogram2D(0)->Draw("COLZ"); m23m13_contour.Draw("P");
	dataDiagrams.getHistogram2D(0)->SetLineWidth(1);
	c1->cd(2);
	fitHitMissDiagrams.getHistogram2D(0)->Draw("COLZ");
	m23m13_contour.Draw("P");
	c1->cd(3);
	fitDiagrams.getHistogram2D(0)->Draw("COLZ");
	m23m13_contour.Draw("P");
	c1->cd(4);
	adaptiveResiduals->Draw("COLZ");
	m23m13_contour.Draw("P");

	//----- plotting invariant mass distributions -----
	TCanvas* c2 = new TCanvas("invmass","invmass",50,50,1200,400);
	c2->Divide(3,1);

	//Plotting mKKsq
	c2->cd(1);
	CreateHist(0);

	//Plotting mKSK+sq
	c2->cd(2);
	CreateHist(1);

	//Plotting mKSK+sq
	c2->cd(3);
	CreateHist(2);
	c2->cd(3);
	TLegend* leg = new TLegend(0.15,0.6,0.50,0.85);
	leg->AddEntry(dataDiagrams.getHistogram(2),"Data");
	leg->AddEntry(fitDiagrams.getHistogram(2),"Model");
	if(ampHistos.size())
		leg->AddEntry(ampHistos.back().getHistogram(2),"Background");
	leg->SetFillStyle(0);
	leg->Draw();

	//----- plotting signal amplitude contributions -----
	TCanvas* c5 = new TCanvas("signalInvmass","signalInvmass",50,50,1200,400);
	c5->Divide(3,1);

	//Plotting mKKsq
	c5->cd(1);
	CreateHist2(0);

	//Plotting mKSK+sq
	c5->cd(2);
	CreateHist2(1);

	//Plotting mKSK+sq
	c5->cd(3);
	CreateHist2(2);
	//	c2->cd(3);
	//	TLegend* leg = new TLegend(0.15,0.6,0.50,0.85);
	//	leg->AddEntry(dataDiagrams.getHistogram(2),"Data");
	//	leg->AddEntry(fitDiagrams.getHistogram(2),"Model");
	//	leg->SetFillStyle(0);
	//	leg->Draw();

	//----- Helicity angle distributions -----
	TCanvas* c3 = new TCanvas("helicityAngle","helicity angle",50,50,1200,800);
	c3->Divide(3,2);
	c3->cd(1); CreateHist(3);
	c3->cd(2); CreateHist(4);
	c3->cd(3); CreateHist(5);
	c3->cd(4); CreateHist(6);
	c3->cd(5); CreateHist(7);
	c3->cd(6); CreateHist(8);

	//----- Weights distributions -----
	TCanvas* c4 = new TCanvas("weights","weights",50,50,1200,400);
	c4->cd(); h_weights.Draw();

	//----- Write to TFile -----
	TFile* tf2 = new TFile(_name+".root","recreate");
	if ( tf2->IsZombie() ) {
		std::cout << "Error opening output file" << std::endl;
		exit(-1);
	}
	m12m13_contour.Write("m12m13_contour",TObject::kOverwrite,0);
	m23m12_contour.Write("m23m12_contour",TObject::kOverwrite,0);
	m23m13_contour.Write("m23m13_contour",TObject::kOverwrite,0);
	c1->Write("dalitz",TObject::kOverwrite,0);
	c2->Write("invmass",TObject::kOverwrite,0);
	c5->Write("signalInvmass",TObject::kOverwrite,0);
	c3->Write("helicityAngle",TObject::kOverwrite,0);

	//Save data trees and histograms
	tf2->mkdir("hist");
	tf2->cd("hist");
	dataAdaptiveBinningHist->Write("dataAdaptive",TObject::kOverwrite,0);
	fitAdaptiveBinningHist->Write("fitAdaptive",TObject::kOverwrite,0);
	adaptiveResiduals->Write("residualsAdaptive",TObject::kOverwrite,0);
	dataDiagrams.Write();
	phspDiagrams.Write();
	fitDiagrams.Write();
	fitHitMissDiagrams.Write();
	auto itAmp = ampHistos.begin();
	for( ; itAmp!=ampHistos.end(); ++itAmp)
		(*itAmp).Write();

	itAmp = signalComponents.begin();
	for( ; itAmp!=signalComponents.end(); ++itAmp)
		(*itAmp).Write();

	//Write some canvas to single files
	c2->Print(_name+"-invmass.root");
	c2->Print(_name+"-invmass.pdf");

	tf2->Close();

	return ;
}

void plotData::CreateHist(unsigned int id)
{
	TPad* pad;
	std::vector<TH1D*> v;
	std::vector<TString> options;
	if( s_data ) {
		v.push_back(dataDiagrams.getHistogram(id));
		options.push_back("E1");
	}
	if( _ampVec.size() ){
		v.push_back(fitDiagrams.getHistogram(id));
		options.push_back("Sames,Hist");
		for(unsigned int i=0; i<plotComponent.size(); ++i){
			ampHistos.at(plotComponent.at(i).first).setColor(plotComponent.at(i).second);
			v.push_back(ampHistos.at(plotComponent.at(i).first).getHistogram(id));
			options.push_back("Sames,Hist");
		}
	}
	pad = drawPull(v, options);
}

void plotData::CreateHist2(unsigned int id)
{
	TPad* pad;
	std::vector<TH1D*> v;
	std::vector<TString> options;
//	v.push_back(ampHistos.at(0).getHistogram(id));
	v.push_back(fitDiagrams.getHistogram(id));
	options.push_back("Hist");

	if( signalComponents.size() ){
		for(unsigned int i=0; i<signalComponents.size(); ++i){
			v.push_back(signalComponents.at(i).getHistogram(id));
			options.push_back("Sames,Hist");
		}
	}
	pad = drawHist(v, options);
}

//===================== dalitzHisto =====================
dalitzHisto::dalitzHisto(std::string n, std::string t,
		unsigned int bins) : name(n), title(t), nBins(bins), _integral(0.0)
{
	ComPWA::Physics::DPKinematics::DalitzKinematics* kin =
			dynamic_cast<ComPWA::Physics::DPKinematics::DalitzKinematics*>(Kinematics::instance());

	//Initialize TTree
	tree = std::unique_ptr<TTree>(
			new TTree( TString(name), TString(title) )
	);

	//Adding branches to TTree
	tree->Branch(TString(name),&t_point);
	tree->Branch("efficiency",&t_eff,"eff/D");
	tree->Branch("weight",&t_weight,"weight/D");

	char label[60];
	for(int i=0; i<kin->GetNVars(); ++i){
		TString varName(kin->GetVarName(i));
		TString varTitle(kin->GetVarTitle(i));

		//Creating TH1D for each variable
		auto limit = kin->GetMinMax(i);

		arr.push_back( TH1D(
				name+varName,
				TString(title+" ")+varTitle,
				nBins,
				limit.first,
				limit.second)
		);
		double binWidth = (double)(limit.second - limit.first)/nBins;
		sprintf(label,"Entries /%f.3",binWidth);
		arr.back().GetYaxis()->SetTitle(label);
		arr.back().GetXaxis()->SetTitle(varTitle);
		arr.back().Sumw2();
	}

	auto m23sq_limit = kin->GetMinMax(0);
	auto m13sq_limit = kin->GetMinMax(1);
	auto m12sq_limit = kin->GetMinMax(2);
	double m23sq_min = m23sq_limit.first;
	double m23sq_max = m23sq_limit.second;
	double m13sq_min = m13sq_limit.first;
	double m13sq_max = m13sq_limit.second;
	double m12sq_min = m12sq_limit.first;
	double m12sq_max = m12sq_limit.second;

	arr2D.push_back( TH2D(TString(name+"_m23sqm13sq"),TString(title),
			nBins,m23sq_min,m23sq_max,nBins,m13sq_min,m13sq_max) );
	arr2D.push_back( TH2D(TString(name+"_m23sqm12sq"),TString(title),
			nBins,m23sq_min,m23sq_max,nBins,m12sq_min,m12sq_max) );
	arr2D.push_back( TH2D(TString(name+"_m12sqm13sq"),TString(title),
			nBins,m12sq_min,m12sq_max,nBins,m13sq_min,m13sq_max) );
	arr2D.push_back( TH2D(TString(name+"_m23sqCosTheta"),TString(title),
			nBins,m23sq_min,m23sq_max,nBins,-1,1) );

	arr2D.at(0).GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}]");
	arr2D.at(0).GetYaxis()->SetTitle("m_{K_{S}K^{+}}^{2} [GeV^{2}/c^{4}]");
	arr2D.at(1).GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}");
	arr2D.at(1).GetYaxis()->SetTitle("m_{K_{S}K^{-}}^{2} [GeV^{2}/c^{4}]");
	arr2D.at(2).GetXaxis()->SetTitle("m_{K_{S}K^{-}}^{2} [GeV^{2}/c^{4}]");
	arr2D.at(2).GetYaxis()->SetTitle("m_{K_{S}K^{+}}^{2} [GeV^{2}/c^{4}]");
	arr2D.at(3).GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}]");
	arr2D.at(3).GetYaxis()->SetTitle("#cos(#Theta)_{KK}");

	auto itr = arr2D.begin();
	for( ; itr!=arr2D.end(); ++itr ){
		(*itr).GetXaxis()->SetNdivisions(508);
		(*itr).GetZaxis()->SetTitle("Entries");
	}
	return;
}

void dalitzHisto::Fill(Event& event, double w)
{
	dataPoint point(event);
	double weight = point.getWeight()*w; //use event weights?

	_integral += weight;
	for(int i=0; i<Kinematics::instance()->GetNVars(); ++i){
		arr.at(i).Fill(point.getVal(i),weight);
	}
	t_eff = event.getEfficiency();
	t_weight = weight;
	t_point = point.getPoint();
	//	tree->Fill();

	double m12sq = point.getVal(2);
	double m13sq = point.getVal(1);
	double m23sq = point.getVal(0);
	double cos12 = point.getVal(4);
	double cos13 = point.getVal(5);
	double cos23 = point.getVal(8);

	arr2D.at(0).Fill(m23sq,m13sq,weight);
	arr2D.at(1).Fill(m23sq,m12sq,weight);
	arr2D.at(2).Fill(m12sq,m13sq,weight);
	arr2D.at(3).Fill(m23sq,cos23,weight);

}

void dalitzHisto::SetStats(bool b)
{
	auto n = arr.size();
	for(int i=0; i<n; ++i){
		arr.at(i).SetStats(b);
	}
	auto n2 = arr2D.size();
	for(int i=0; i<n2; ++i){
		arr2D.at(i).SetStats(b);
	}
}

void dalitzHisto::Scale(double w)
{
	auto n = arr.size();
	for(int i=0; i<n; ++i){
		arr.at(i).Scale(w);
	}
	auto n2 = arr2D.size();
	for(int i=0; i<n2; ++i){
		arr2D.at(i).Scale(w);
	}
}

void dalitzHisto::setColor(Color_t color)
{
	auto n = arr.size();
	for(int i=0; i<n; ++i){
		arr.at(i).SetLineColor(color);
		arr.at(i).SetMarkerColor(color);
	}
}

TH1D* dalitzHisto::getHistogram(unsigned int num)
{
	return &arr.at(num);
}

TH2D* dalitzHisto::getHistogram2D(unsigned int num)
{
	return &arr2D.at(num);
}

TH2Poly* dalitzHisto::getTH2PolyPull(TH2Poly* hist1, TH2Poly* hist2, TString name)
{
	if(hist1->GetBins()->GetEntries() != hist2->GetBins()->GetEntries() ){
		std::cout<<"binning doesnt match"<<std::endl;
		return 0;
	}
	TH2Poly* resHist = (TH2Poly*)hist1->Clone(name+hist1->GetName());
	resHist->Reset("");

	hist2->Scale(hist1->Integral()/hist2->Integral());
	double limit=0;
	for(int bin=1; bin<=hist1->GetBins()->GetEntries(); bin++) {
		double c1 = hist1->GetBinContent(bin);
		double c2 = hist2->GetBinContent(bin);
		double c1Err = hist1->GetBinError(bin);
		double c2Err = hist2->GetBinError(bin);
		if(c1>0 && c2>0) {
			double pull =(c1-c2)/sqrt(c1Err*c1Err+c2Err*c2Err);
			resHist->SetBinContent(bin,pull);
			if(std::abs(pull)>limit) limit=std::abs(pull);
		}
		else
			resHist->SetBinContent(bin,-999);//empty bins are set to error value
		resHist->SetBinError(bin,0);
	}
	limit=4;
	limit+=1;
	resHist->GetZaxis()->SetRangeUser(-limit,limit);//symmetric range including all values
	resHist->GetZaxis()->SetTitle("deviation [#sigma]");
	resHist->GetZaxis()->SetTitleOffset(0.5);

	return resHist;
}

void dalitzHisto::Write()
{
	tree->Write(TString(name)+"_tree");
	gDirectory->mkdir(TString(name)+"_hist");
	gDirectory->cd(TString(name)+"_hist");
	auto n = arr.size();
	for(int i=0; i<n; ++i){
		arr.at(i).Write();
	}
	auto n2 = arr2D.size();
	for(int i=0; i<n2; ++i){
		arr2D.at(i).Write();
	}
	gDirectory->cd("../");
}
