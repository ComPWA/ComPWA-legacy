#include "Analysis/PlotData.hpp"

plotData::plotData(std::string name, std::string file,std::shared_ptr<Data> set1,std::shared_ptr<Data> set2,std::shared_ptr<Amplitude> amp):
_set1(set1),_set2(set2),_amp(amp){
	_outFile=file;
	_name=name;
	_style=0;
};
plotData::plotData(std::string name, std::string file,std::shared_ptr<Data> set1,std::shared_ptr<Data> set2):
				_set1(set1),_set2(set2),_par(){
	_outFile=file;
	_name=name;
	_style=0;
};
plotData::plotData(std::string name, std::string file,std::shared_ptr<Data> set1):
				_set1(set1),_set2(std::shared_ptr<Data>(new RootReader())),_par(){
	_outFile=file;
	_name=name;
	_style=1;
};
void plotData::setStyle(unsigned int style){
	_style=style;
};

void plotData::plot(){

	double min_range=.9;
	double max_range=2.;
	int nBins=100;
	//Plot result
	TH2D* bw2312 = new TH2D("bw2312","DATA",nBins,min_range,max_range,nBins,min_range,max_range);
	bw2312->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}");
	bw2312->GetXaxis()->CenterTitle();
	bw2312->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
	bw2312->GetYaxis()->CenterTitle();
	bw2312->SetLineColor(kBlack);
	TH2D* bw1213 = new TH2D("bw1213","DATA",nBins,min_range,max_range,nBins,min_range,max_range);
	bw1213->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
	bw1213->GetXaxis()->CenterTitle();
	bw1213->GetYaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
	bw1213->GetYaxis()->CenterTitle();
	bw1213->SetLineColor(kBlack);
	TH2D* bw2313 = new TH2D("bw2313","DATA",nBins,min_range,max_range,nBins,min_range,max_range);
	bw2313->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}");
	bw2313->GetXaxis()->CenterTitle();
	bw2313->GetYaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
	bw2313->GetYaxis()->CenterTitle();
	bw2313->SetLineColor(kBlack);

	//double masssq12PHSP, masssq13PHSP, masssq23PHSP;
	TH2D* bw2312PHSP = new TH2D("bw2312PHSP","PHSP",nBins,min_range,max_range,nBins,min_range,max_range);
	bw2312PHSP->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}");
	bw2312PHSP->GetXaxis()->CenterTitle();
	bw2312PHSP->GetYaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
	bw2312PHSP->GetYaxis()->CenterTitle();
	bw2312PHSP->SetLineColor(kRed);
	TH2D* bw1213PHSP = new TH2D("bw1213PHSP","PHSP",nBins,min_range,max_range,nBins,min_range,max_range);
	bw1213PHSP->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
	bw1213PHSP->GetXaxis()->CenterTitle();
	bw1213PHSP->GetYaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
	bw1213PHSP->GetYaxis()->CenterTitle();
	bw1213PHSP->SetLineColor(kRed);
	TH2D* bw2313PHSP = new TH2D("bw2313PHSP","PHSP",nBins,min_range,max_range,nBins,min_range,max_range);
	bw2313PHSP->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}");
	bw2313PHSP->GetXaxis()->CenterTitle();
	bw2313PHSP->GetYaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
	bw2313PHSP->GetYaxis()->CenterTitle();
	bw2313PHSP->SetLineColor(kRed);

	TH2D* bw2312FIT = new TH2D("bw2312FIT","FIT",nBins,min_range,max_range,nBins,min_range,max_range);
	bw2312FIT->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}");
	bw2312FIT->GetXaxis()->CenterTitle();
	bw2312FIT->GetYaxis()->SetTitle("#m_{12}^{2} / GeV^{2}");
	bw2312FIT->GetYaxis()->CenterTitle();
	bw2312FIT->SetLineColor(kBlue);
	TH2D* bw1213FIT = new TH2D("bw1213FIT","FIT",nBins,min_range,max_range,nBins,min_range,max_range);
	bw1213FIT->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
	bw1213FIT->GetXaxis()->CenterTitle();
	bw1213FIT->GetYaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
	bw1213FIT->GetYaxis()->CenterTitle();
	bw1213FIT->SetLineColor(kBlue);
	TH2D* bw2313FIT = new TH2D("bw2313FIT","FIT",nBins,min_range,max_range,nBins,min_range,max_range);
	bw2313FIT->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}");
	bw2313FIT->GetXaxis()->CenterTitle();
	bw2313FIT->GetYaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
	bw2313FIT->GetYaxis()->CenterTitle();
	bw2313FIT->SetLineColor(kBlue);

	double masssq12, masssq13, masssq23;
	for(unsigned int i = 0; i < _set1->getNEvents(); i++){
		Event event(_set1->getEvent(i));

		//myReader.getEvent(-1, a, b, masssq);
		//if(!myReader.getEvent(i, event)) continue; TODO: try exception
		if(!event.getNParticles() == 3) continue;
		//if(!event) continue;
		//cout << "Event: \t" << i << "\t NParticles: \t" << event.getNParticles() << endl;
		const Particle &a(event.getParticle(0));
		const Particle &b(event.getParticle(1));
		const Particle &c(event.getParticle(2));
		masssq12 = pow(a.E+b.E,2) - pow(a.px+b.px ,2) - pow(a.py+b.py ,2) - pow(a.pz+b.pz ,2);
		masssq13 = pow(a.E+c.E,2) - pow(a.px+c.px ,2) - pow(a.py+c.py ,2) - pow(a.pz+c.pz ,2);
		masssq23 = pow(b.E+c.E,2) - pow(b.px+c.px ,2) - pow(b.py+c.py ,2) - pow(b.pz+c.pz ,2);

		bw2312->Fill(masssq23,masssq12);
		bw1213->Fill(masssq12,masssq13);
		bw2313->Fill(masssq23,masssq13);
	}

	for(unsigned int i = 0; i < _set2->getNEvents(); i++){
		Event event(_set2->getEvent(i));

		//myReader.getEvent(-1, a, b, masssq);
		//if(!myReader.getEvent(i, event)) continue; TODO: try exception
		if(!event.getNParticles() == 3) continue;
		//if(!event) continue;
		//cout << "Event: \t" << i << "\t NParticles: \t" << event.getNParticles() << endl;
		const Particle &a(event.getParticle(0));
		const Particle &b(event.getParticle(1));
		const Particle &c(event.getParticle(2));
		masssq12 = pow(a.E+b.E,2) - pow(a.px+b.px ,2) - pow(a.py+b.py ,2) - pow(a.pz+b.pz ,2);
		masssq13 = pow(a.E+c.E,2) - pow(a.px+c.px ,2) - pow(a.py+c.py ,2) - pow(a.pz+c.pz ,2);
		masssq23 = pow(b.E+c.E,2) - pow(b.px+c.px ,2) - pow(b.py+c.py ,2) - pow(b.pz+c.pz ,2);

		bw2312PHSP->Fill(masssq23,masssq12);
		bw1213PHSP->Fill(masssq12,masssq13);
		bw2313PHSP->Fill(masssq23,masssq13);

		vector<double> x;
		x.push_back(sqrt(masssq23));
		x.push_back(sqrt(masssq13));
//		x.push_back(sqrt(masssq12));
		if(_par.GetNParameter()!=0){
			bw2312FIT->Fill(masssq23,masssq12,_amp->intensity(x,_par));
			bw1213FIT->Fill(masssq12,masssq13,_amp->intensity(x,_par));
			bw2313FIT->Fill(masssq23,masssq13,_amp->intensity(x,_par));
		}
	}
	TCanvas* c1 = new TCanvas(_name.c_str(),_name.c_str(),200,10,1400,700);
	switch(_style){
	case 0:
		c1->Divide(3,2);
		c1->cd(1); bw2313->Draw("COLZ");
		c1->cd(2); bw2313FIT->Draw("COLZ");
		c1->cd(3); bw2313PHSP->Draw("COLZ");
		c1->cd(4); bw2313->ProjectionX()->Draw();bw2313FIT->ProjectionX()->Draw("Sames");bw2313PHSP->ProjectionX()->Draw("Sames");
		gPad->BuildLegend();
		c1->cd(5); bw2313->ProjectionY()->Draw();bw2313FIT->ProjectionY()->Draw("Sames");bw2313PHSP->ProjectionY()->Draw("Sames");
		gPad->BuildLegend();
		c1->cd(6); bw1213->ProjectionX()->Draw();bw1213FIT->ProjectionX()->Draw("Sames");bw1213PHSP->ProjectionX()->Draw("Sames");
		gPad->BuildLegend();
		break;
	case 1://style for one dataset only
		c1->Divide(3,2);
		c1->cd(1); bw2313->Draw("COLZ");
		c1->cd(2); bw2312->Draw("COLZ");
		c1->cd(3); bw1213->Draw("COLZ");
		c1->cd(4); bw2313->ProjectionX()->Draw();
		c1->cd(5); bw2313->ProjectionY()->Draw();
		c1->cd(6); bw1213->ProjectionX()->Draw();
		break;
	default:
		cout<<"Choose a style for plotting!"<<endl;
	}

	TFile* tf = new TFile(_outFile.c_str(),"update");
	if ( tf->IsZombie() ) {
		std::cout << "Error opening output file" << std::endl;
		exit(-1);
	}
	c1->Write("",TObject::kOverwrite,0);
	tf->Write();
	tf->Close();
	bw2312->Delete();
	bw1213->Delete();
	bw2313->Delete();
	bw2312PHSP->Delete();
	bw1213PHSP->Delete();
	bw2313PHSP->Delete();
	bw2312FIT->Delete();
	bw1213FIT->Delete();
	bw2313FIT->Delete();
	//	c1->Delete();
	return ;

};
