#ifndef HISTTOOLS_H
#define HISTTOOLS_H 

#include <TGraph.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TPad.h>
#include <TMath.h>
#include <TKDTreeBinning.h>

Double_t Chi2TestX(const TH1* h1, const TH1* h2,  Double_t &chi2, Int_t &ndf,
		Int_t &igood, Option_t *option="",  Double_t *res=0)
{
	// The computation routine of the Chisquare test. For the method description,
	// see Chi2Test() function.
	// Returns p-value
	// parameters:
	//  - h2-second histogram
	//  - option:
	//     "UU" = experiment experiment comparison (unweighted-unweighted)
	//     "UW" = experiment MC comparison (unweighted-weighted). Note that the first
	//           histogram should be unweighted
	//     "WW" = MC MC comparison (weighted-weighted)
	//
	//     "NORM" = if one or both histograms is scaled
	//
	//     "OF" = overflows included
	//     "UF" = underflows included
	//         by default underflows and overflows are not included
	//
	//  - igood:
	//       igood=0 - no problems
	//        For unweighted unweighted  comparison
	//       igood=1'There is a bin in the 1st histogram with less than 1 event'
	//       igood=2'There is a bin in the 2nd histogram with less than 1 event'
	//       igood=3'when the conditions for igood=1 and igood=2 are satisfied'
	//        For  unweighted weighted  comparison
	//       igood=1'There is a bin in the 1st histogram with less then 1 event'
	//       igood=2'There is a bin in the 2nd histogram with less then 10 effective number of events'
	//       igood=3'when the conditions for igood=1 and igood=2 are satisfied'
	//        For  weighted weighted  comparison
	//       igood=1'There is a bin in the 1st  histogram with less then 10 effective
	//        number of events'
	//       igood=2'There is a bin in the 2nd  histogram with less then 10 effective
	//               number of events'
	//       igood=3'when the conditions for igood=1 and igood=2 are satisfied'
	//
	//  - chi2 - chisquare of the test
	//  - ndf  - number of degrees of freedom (important, when both histograms have the same
	//         empty bins)
	//  - res -  normalized residuals for further analysis


	Int_t i = 0, j=0, k = 0;
	Int_t i_start, i_end;
	Int_t j_start, j_end;
	Int_t k_start, k_end;

	Double_t bin1, bin2;
	Double_t err1,err2;
	Double_t sum1=0, sum2=0;
	Double_t sumw1=0, sumw2=0;


	chi2 = 0;
	ndf = 0;

	TString opt = option;
	opt.ToUpper();

	const TAxis *xaxis1 = h1->GetXaxis();
	const TAxis *xaxis2 = h2->GetXaxis();
	const TAxis *yaxis1 = h1->GetYaxis();
	const TAxis *yaxis2 = h2->GetYaxis();
	const TAxis *zaxis1 = h1->GetZaxis();
	const TAxis *zaxis2 = h2->GetZaxis();

	Int_t nbinx1 = xaxis1->GetNbins();
	Int_t nbinx2 = xaxis2->GetNbins();
	Int_t nbiny1 = yaxis1->GetNbins();
	Int_t nbiny2 = yaxis2->GetNbins();
	Int_t nbinz1 = zaxis1->GetNbins();
	Int_t nbinz2 = zaxis2->GetNbins();

	//check dimensions
	if (h1->GetDimension() != h2->GetDimension() ){
		Error("Chi2TestX","Histograms have different dimensions.");
		return 0;
	}

	//check number of channels
	if (nbinx1 != nbinx2) {
		Error("Chi2TestX","different number of x channels");
	}
	if (nbiny1 != nbiny2) {
		Error("Chi2TestX","different number of y channels");
	}
	if (nbinz1 != nbinz2) {
		Error("Chi2TestX","different number of z channels");
	}

	//check for ranges
	i_start = j_start = k_start = 1;
	i_end = nbinx1;
	j_end = nbiny1;
	k_end = nbinz1;

	if (xaxis1->TestBit(TAxis::kAxisRange)) {
		i_start = xaxis1->GetFirst();
		i_end   = xaxis1->GetLast();
	}
	if (yaxis1->TestBit(TAxis::kAxisRange)) {
		j_start = yaxis1->GetFirst();
		j_end   = yaxis1->GetLast();
	}
	if (zaxis1->TestBit(TAxis::kAxisRange)) {
		k_start = zaxis1->GetFirst();
		k_end   = zaxis1->GetLast();
	}


	if (opt.Contains("OF")) {
		if (h1->GetDimension() == 3) k_end = ++nbinz1;
		if (h1->GetDimension() >= 2) j_end = ++nbiny1;
		if (h1->GetDimension() >= 1) i_end = ++nbinx1;
	}

	if (opt.Contains("UF")) {
		if (h1->GetDimension() == 3) k_start = 0;
		if (h1->GetDimension() >= 2) j_start = 0;
		if (h1->GetDimension() >= 1) i_start = 0;
	}

	ndf = (i_end - i_start + 1)*(j_end - j_start + 1)*(k_end - k_start + 1) - 1;

	Bool_t comparisonUU = opt.Contains("UU");
	Bool_t comparisonUW = opt.Contains("UW");
	Bool_t comparisonWW = opt.Contains("WW");
	Bool_t scaledHistogram  = opt.Contains("NORM");
	if (scaledHistogram && !comparisonUU) {
		Info("Chi2TestX","NORM option should be used together with UU option. It is ignored");
	}
	// look at histo global bin content and effective entries
	const Int_t kNstat=30;
	Stat_t s[kNstat];
	h1->GetStats(s);// s[1] sum of squares of weights, s[0] sum of weights
	double sumBinContent1 = s[0];
	double effEntries1 = (s[1] ? s[0]*s[0]/s[1] : 0.);

	h2->GetStats(s);// s[1] sum of squares of weights, s[0] sum of weights
	double sumBinContent2 = s[0];
	double effEntries2 = (s[1] ? s[0]*s[0]/s[1] : 0.);

	if (!comparisonUU && !comparisonUW && !comparisonWW ) {
		// deduce automatically from type of histogram
		if (TMath::Abs(sumBinContent1 - effEntries1) < 1) {
			if ( TMath::Abs(sumBinContent2 - effEntries2) < 1) comparisonUU = true;
			else comparisonUW = true;
		}
		else comparisonWW = true;
	}
	// check unweighted histogram
	if (comparisonUW) {
		if (TMath::Abs(sumBinContent1 - effEntries1) >= 1) {
			Warning("Chi2TestX","First histogram is not unweighted and option UW has been requested");
		}
	}
	if ( (!scaledHistogram && comparisonUU)   ) {
		if ( ( TMath::Abs(sumBinContent1 - effEntries1) >= 1) || (TMath::Abs(sumBinContent2 - effEntries2) >= 1) ) {
			Warning("Chi2TestX","Both histograms are not unweighted and option UU has been requested");
		}
	}


	//get number of events in histogram
	if (comparisonUU && scaledHistogram) {
		for (i=i_start; i<=i_end; i++) {
			for (j=j_start; j<=j_end; j++) {
				for (k=k_start; k<=k_end; k++) {
					Int_t bin = h1->GetBin(i,j,k);
					bin1 = h1->GetBinContent(bin);
					bin2 = h2->GetBinContent(bin);
					err1 = h1->GetBinError(bin);
					err2 = h2->GetBinError(bin);
					if (err1 > 0 ) {
						bin1 *= bin1/(err1*err1);
						//avoid rounding errors
						bin1 = TMath::Floor(bin1+0.5);
					}
					else
						bin1 = 0;

					if (err2 > 0) {
						bin2 *= bin2/(err2*err2);
						//avoid rounding errors
						bin2 = TMath::Floor(bin2+0.5);
					}
					else
						bin2 = 0;

					// sum contents
					sum1 += bin1;
					sum2 += bin2;
					sumw1 += err1*err1;
					sumw2 += err2*err2;
				}
			}
		}
		if (sumw1 <= 0 || sumw2 <= 0) {
			Error("Chi2TestX","Cannot use option NORM when one histogram has all zero errors");
			return 0;
		}

	} else {
		for (i=i_start; i<=i_end; i++) {
			for (j=j_start; j<=j_end; j++) {
				for (k=k_start; k<=k_end; k++) {
					Int_t bin = h1->GetBin(i,j,k);
					//			   if(j==0)
					// std::cout<<i<<" "<<j<<" "<<k<< " -> "<<bin<<" "
					// <<h1->GetBinContent(bin)<<" "
					// <<h2->GetBinContent(bin)<< std::endl;
					sum1 += h1->GetBinContent(bin);
					sum2 += h2->GetBinContent(bin);
					if ( comparisonWW ) {
						err1 = h1->GetBinError(bin);
						sumw1 += err1*err1;
					}
					if ( comparisonUW || comparisonWW ) {
						err2 = h2->GetBinError(bin);
						sumw2 += err2*err2;
					}
				}
			}
		}
	}
	//  std::cout<<sum1<<" "<<sum2<<std::endl;
	//checks that the histograms are not empty
	if (sum1 == 0 || sum2 == 0) {
		Error("Chi2TestX","one histogram is empty");
		return 0;
	}

	if ( comparisonWW  && ( sumw1 <= 0 && sumw2 <=0 ) ){
		Error("Chi2TestX","Hist1 and Hist2 have both all zero errors\n");
		return 0;
	}

	//THE TEST
	Int_t m=0, n=0;

	//Experiment - experiment comparison
	if (comparisonUU) {
		Double_t sum = sum1 + sum2;
		for (i=i_start; i<=i_end; i++) {
			for (j=j_start; j<=j_end; j++) {
				for (k=k_start; k<=k_end; k++) {
					Int_t bin = h1->GetBin(i,j,k);
					bin1 = h1->GetBinContent(bin);
					bin2 = h2->GetBinContent(bin);


					if (scaledHistogram) {
						// scale bin value to effective bin entries
						err1 = h1->GetBinError(bin);
						if (err1 > 0 ) {
							bin1 *= bin1/(err1*err1);
							//avoid rounding errors
							bin1 = TMath::Floor(bin1+0.5);
						}
						else
							bin1 = 0;

						err2 = h2->GetBinError(bin);
						if (err2 > 0) {
							bin2 *= bin2/(err2*err2);
							//avoid rounding errors
							bin2 = TMath::Floor(bin2+0.5);
						}
						else
							bin2 = 0;

					}

					if ( (int(bin1) == 0)  && (int(bin2) == 0) ) {
						--ndf;  //no data means one degree of freedom less
					} else {


						Double_t binsum = bin1 + bin2;
						Double_t nexp1 = binsum*sum1/sum;
						//Double_t nexp2 = binsum*sum2/sum;

						//                  if(opt.Contains("P")) printf("bin %d p = %g\t",i,binsum/sum);

						if (res)
							res[i-i_start] = (bin1-nexp1)/TMath::Sqrt(nexp1);

						if (bin1 < 1) {
							m++;
						}
						if (bin2 < 1) {
							n++;
						}

						//Habermann correction for residuals
						Double_t correc = (1-sum1/sum)*(1-binsum/sum);
						if (res) {
							res[i-i_start] /= TMath::Sqrt(correc);
						}

						Double_t delta = sum2*bin1-sum1*bin2;
						chi2 += delta*delta/binsum;

					}
				}
			}
		}

		chi2 /= (sum1*sum2);
		// flag error only when of the two histogram is zero
		if (m) {
			igood += 1;
			Info("Chi2TestX","There is a bin in h1 with less than 1 event.\n");
		}
		if (n) {
			igood += 2;
			Info("Chi2TestX","There is a bin in h2 with less than 1 event.\n");
		}

		Double_t prob = TMath::Prob(chi2,ndf);
		return prob;

	}


	//unweighted - weighted  comparison
	// case of err2 = 0 and bin2 not zero is treated without problems
	// by excluding second chi2 sum
	// and can be considered as a comparison data-theory
	if ( comparisonUW ) {
		for (i=i_start; i<=i_end; i++) {
			for (j=j_start; j<=j_end; j++) {
				for (k=k_start; k<=k_end; k++) {
					Int_t x=0;
					Int_t bin = h1->GetBin(i,j,k);
					bin1 = h1->GetBinContent(bin);
					bin2 = h2->GetBinContent(bin);
					err2 = h2->GetBinError(bin);

					err2 *= err2;

					// case both histogram have zero bin contents
					if ( (int(bin1) == 0) && (bin2*bin2 == 0) ) {
						--ndf;  //no data means one degree of freedom less
						continue;
					}

					// case weighted histogram has zero bin content and error
					if (bin2*bin2 == 0 && err2 == 0) {
						if (sumw2 > 0) {
							// use as approximated  error as 1 scaled by a scaling ratio
							// estimated from the total sum weight and sum weight squared
							err2 = sumw2/sum2;
						}
						else {
							// return error because infinite discrepancy here:
							// bin1 != 0 and bin2 =0 in a histogram with all errors zero
							Error("Chi2TestX","Hist2 has in bin %d,%d,%d zero content and all zero errors\n", i,j,k);
							chi2 = 0; return 0;
						}
					}

					if (bin1 < 1)  m++;
					if (err2 > 0 && bin2*bin2/err2 < 10) n++;

					Double_t var1 = sum2*bin2 - sum1*err2;
					Double_t var2 = var1*var1 + 4*sum2*sum2*bin1*err2;
					// if bin1 is zero and bin2=1 and sum1=sum2 var1=0 && var2 ==0
					// approximate by adding +1 to bin1
					// LM (h1 need to be fixed for numerical errors)
					while (var1*var1+bin1 == 0 || var1+var2 == 0) {
						sum1++;
						bin1++;
						x++;
						var1 = sum2*bin2 - sum1*err2;
						var2 = var1*var1 + 4*sum2*sum2*bin1*err2;
					}
					var2 = TMath::Sqrt(var2);
					while (var1+var2 == 0) {
						sum1++;
						bin1++;
						x++;
						var1 = sum2*bin2 - sum1*err2;
						var2 = var1*var1 + 4*sum2*sum2*bin1*err2;
						while (var1*var1+bin1 == 0 || var1+var2 == 0) {
							sum1++;
							bin1++;
							x++;
							var1 = sum2*bin2 - sum1*err2;
							var2 = var1*var1 + 4*sum2*sum2*bin1*err2;
						}
						var2 = TMath::Sqrt(var2);
					}

					Double_t probb = (var1+var2)/(2*sum2*sum2);

					Double_t nexp1 = probb * sum1;
					Double_t nexp2 = probb * sum2;

					//               if(opt.Contains("P")) printf("bin %d p = %g\t",i,probb);

					Double_t delta1 = bin1 - nexp1;
					Double_t delta2 = bin2 - nexp2;

					chi2 += delta1*delta1/nexp1;


					if (err2 > 0) {
						//            	   std::cout<<(bin1-bin2)*(bin1-bin2)/nexp1
						//            			   <<delta1*delta1/nexp1+delta2*delta2/err2
						//            			   <<std::endl;
						chi2 += delta2*delta2/err2;
					}

					if (res) {
						if (err2 > 0) {
							Double_t temp1 = sum2*err2/var2;
							Double_t temp2 = 1 + (sum1*err2 - sum2*bin2)/var2;
							temp2 = temp1*temp1*sum1*probb*(1-probb) + temp2*temp2*err2/4;
							// invert sign here
							res[i-i_start] = - delta2/TMath::Sqrt(temp2);
						}
						else
							res[i-i_start] = delta1/TMath::Sqrt(nexp1);

					}
				}
			}
		}

		if (m) {
			igood += 1;
			Info("Chi2TestX","There is a bin in h1 with less than 1 event.\n");
		}
		if (n) {
			igood += 2;
			Info("Chi2TestX","There is a bin in h2 with less than 10 effective events.\n");
		}

		Double_t prob = TMath::Prob(chi2,ndf);

		return prob;
	}

	// weighted - weighted  comparison
	if (comparisonWW) {
		for (i=i_start; i<=i_end; i++) {
			for (j=j_start; j<=j_end; j++) {
				for (k=k_start; k<=k_end; k++) {
					Int_t bin = h1->GetBin(i,j,k);
					bin1 = h1->GetBinContent(bin);
					bin2 = h2->GetBinContent(bin);
					err1 = h1->GetBinError(bin);
					err2 = h2->GetBinError(bin);
					err1 *= err1;
					err2 *= err2;

					// case both histogram have zero bin contents
					// (use square of bin1 to avoid numerical errors)
					if ( (bin1*bin1 == 0) && (bin2*bin2 == 0) ) {
						--ndf;  //no data means one degree of freedom less
						continue;
					}

					if ( (err1 == 0) && (err2 == 0) ) {
						// cannot treat case of booth histogram have zero zero errors
						Error("Chi2TestX","h1 and h2 both have bin %d,%d,%d with all zero errors\n", i,j,k);
						chi2 = 0; return 0;
					}

					Double_t sigma  = sum1*sum1*err2 + sum2*sum2*err1;
					Double_t delta = sum2*bin1 - sum1*bin2;
					chi2 += delta*delta/sigma;

					//               if(opt.Contains("P")) printf("bin %d p = %g\t",i, (bin1*sum1/err1 + bin2*sum2/err2)/(sum1*sum1/err1 + sum2*sum2/err2));

					if (res) {
						Double_t temp = bin1*sum1*err2 + bin2*sum2*err1;
						Double_t probb = temp/sigma;
						Double_t z = 0;
						if (err1 > err2 ) {
							Double_t d1 = (bin1 - sum1 * probb);
							Double_t s1 = err1* ( 1. - err2 * sum1 * sum1 / sigma );
							z = d1/ TMath::Sqrt(s1);
						}
						else {
							Double_t d2 = (bin2 - sum2 * probb);
							Double_t s2 = err2* ( 1. - err1 * sum2 * sum2 / sigma );
							z = -d2/ TMath::Sqrt(s2);
						}

						res[i-i_start] = z;
					}

					if (err1 > 0 && bin1*bin1/err1 < 10) m++;
					if (err2 > 0 && bin2*bin2/err2 < 10) n++;
				}
			}
		}
		if (m) {
			igood += 1;
			Info("Chi2TestX","There is a bin in h1 with less than 10 effective events.\n");
		}
		if (n) {
			igood += 2;
			Info("Chi2TestX","There is a bin in h2 with less than 10 effective events.\n");
		}
		Double_t prob = TMath::Prob(chi2,ndf);
		return prob;
	}
	return 0;
}

TH1* getPull(TH1* hist1, TH1* hist2, TString name="pull_")
{
	if(hist1->GetNbinsX() != hist2->GetNbinsX() ||
			hist1->GetNbinsY() != hist2->GetNbinsY()||
			hist1->GetNbinsZ() != hist2->GetNbinsZ()) {
		std::cout<<"binning doesnt match"<<std::endl;
		return 0;
	}
	TH1* pullHist = (TH1*)hist1->Clone(name+hist1->GetName());
	pullHist->Reset();
	double limitPull=0;
	for(int i=1; i<=hist1->GetNbinsX(); i++) {
		for(int j=1; j<=hist1->GetNbinsY(); j++) {
			for(int k=1; k<=hist1->GetNbinsZ(); k++) {
				unsigned int bin = hist1->GetBin(i,j,k);
				double c1 = hist1->GetBinContent(bin);
				double c2 = hist2->GetBinContent(bin);
				double c1Err = hist1->GetBinError(bin);
				double c2Err = hist2->GetBinError(bin);
				double pull;
				if(c1>0 && c2>0) {
					pull =(c1-c2)/sqrt(c1Err*c1Err+c2Err*c2Err);
					pullHist->SetBinContent(bin,pull);
					if(std::fabs(pull)>limitPull) limitPull=abs(pull);
				}
				else
					pullHist->SetBinContent(bin,-999);//empty bins are set to error value
				pullHist->SetBinError(bin,.0001);
			}
		}
	}
	limitPull=4; //set fix limits
	limitPull+=1;
	//X
	pullHist->GetXaxis()->SetTitleSize(.14);
	pullHist->GetXaxis()->SetTitleOffset(.93);
	pullHist->GetXaxis()->SetLabelSize(.12);
	if(pullHist->GetDimension()==1){
		pullHist->GetYaxis()->SetRangeUser(-limitPull,limitPull);
		pullHist->GetYaxis()->SetTitle("deviation [#sigma]");
		pullHist->GetYaxis()->SetTitleOffset(0.36);
		pullHist->GetYaxis()->SetTitleSize(0.12);
		pullHist->GetYaxis()->SetLabelSize(0.12);
		pullHist->GetYaxis()->SetNdivisions(504);
		pullHist->GetYaxis()->CenterTitle();
	}
	if(pullHist->GetDimension()==2){//symmetric range including all values
		pullHist->GetZaxis()->SetRangeUser(-limitPull,limitPull);
		pullHist->GetZaxis()->SetTitle("deviation [#sigma]");
		pullHist->GetYaxis()->SetTitleSize(.14);
		pullHist->GetYaxis()->SetTitleOffset(.93);
		pullHist->GetYaxis()->SetLabelSize(.12);
	}
	pullHist->SetTitle("");
	pullHist->SetStats(0);
	return pullHist;
}
void ttttt()
{
	std::cout<<"1231231231"<<std::endl;
}

TH1* getResidual(TH1* hist1, TH1* hist2, TString name="res_")
{
	if(hist1->GetNbinsX() != hist2->GetNbinsX() ||
			hist1->GetNbinsY() != hist2->GetNbinsY()||
			hist1->GetNbinsZ() != hist2->GetNbinsZ()) {
		std::cout<<"binning doesnt match"<<std::endl;
		return 0;
	}
	TH1* resHist = (TH1*)hist1->Clone(name+hist1->GetName());
	resHist->Reset();

	hist2->Scale(hist1->Integral()/hist2->Integral());
	double limit=0;
	for(int i=1; i<=hist1->GetNbinsX(); i++) {
		for(int j=1; j<=hist1->GetNbinsY(); j++) {
			for(int k=1; k<=hist1->GetNbinsZ(); k++) {
				unsigned int bin = hist1->GetBin(i,j,k);
				double c1 = hist1->GetBinContent(bin);
				double c2 = hist2->GetBinContent(bin);
				if( c1==0 && c2==0 ) {
					double res = -9999;
					resHist->SetBinContent(bin,res);
				} else {
					double res = c1-c2;
					resHist->SetBinContent(bin,res);
					if(std::fabs(res)>limit) limit=abs(res);
				}
				resHist->SetBinError(bin,0);
			}
		}
	}
	limit+=1;
	if(resHist->GetDimension()==1)
		resHist->GetYaxis()->SetRangeUser(-limit,limit);//symmetric range including all values
	if(resHist->GetDimension()==2)
		resHist->GetZaxis()->SetRangeUser(-limit,limit);//symmetric range including all values
	return resHist;
}

TPad* drawHist(std::vector<TH1D*> hist, std::vector<TString> drawOption,
		double min=0, double max=0)
{
	if(hist.size() != drawOption.size() )
		throw std::runtime_error("drawPull() | Number of histograms and number "
				"of draw options does not match!");

	Int_t optTitle = gStyle->GetOptTitle();
	Int_t optStat = gStyle->GetOptStat();
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);//entries only

	TPad* pad = new TPad();
	if( !hist.size() ) return 0;
	for(int i=0; i<hist.size(); i++)
		hist.at(i)->Draw(drawOption.at(i));
	return pad;
}

TPad* drawPull(std::vector<TH1D*> hist, std::vector<TString> drawOption,
		double min=0, double max=0)
{
	if(hist.size() != drawOption.size() )
		throw std::runtime_error("drawPull() | Number of histograms and number "
				"of draw options does not match!");

	Int_t optTitle = gStyle->GetOptTitle();
	Int_t optStat = gStyle->GetOptStat();
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);//entries only

	TPad* pad;
	if( !hist.size() ) return 0;
	else if( hist.size() == 1 ) {
		TPad* pad = new TPad();
		hist.at(0)->Draw(drawOption.at(0));
	} else if(hist.size() > 1){
		TPad* pad = new TPad("hist","hist",0.0,0.3,1,1);
		pad->Draw(); pad->SetMargin(0.1,0.05,0.0,0.05);
		TPad* pad_pull = new TPad("pull","pull",0.0,0.0,1,0.3);
		pad_pull->Draw(); pad_pull->SetMargin(0.1,0.05,0.3,0.0);//left-right-bottom-top

		if( hist.at(0)->GetDimension() != 1 || hist.at(1)->GetDimension() != 1 ) {
			std::cout<<"Dimension of histograms larger 1!"<<std::endl;
			return pad;
		}
		if( hist.at(0)->GetNbinsX() != hist.at(1)->GetNbinsX() ) {
			std::cout<<"binning doesnt match"<<std::endl;
			return pad;
		}

		pad->cd();
		hist.at(0)->GetYaxis()->SetRangeUser(0.00000001, hist.at(0)->GetMaximum()*1.3);
		hist.at(0)->GetYaxis()->SetTitleOffset(0.83);

		hist.at(0)->Draw(drawOption.at(0));
		hist.at(1)->Draw(drawOption.at(1));
		for(int i=2; i<hist.size(); ++i)
			hist.at(i)->Draw(drawOption.at(i));

		Double_t chi2;
		Double_t res[hist.at(0)->GetNbinsX()];
		Int_t ndf, igood;
		hist.at(0)->Chi2TestX(hist.at(1),chi2,ndf,igood,"UW",res);

		char chi2Char[60];
		sprintf(chi2Char,"#chi^{2}_{1D}/ndf = %.2f/%d",chi2,ndf);
		TLatex* ltx = new TLatex();
		ltx->DrawLatexNDC(0.2,0.85,chi2Char);//needs to be drawn inside a TPad
		delete ltx;

		TAxis *xaxis = ((TH1*)hist.at(0))->GetXaxis();
		Int_t fNpoints = xaxis->GetNbins();
		Double_t fX[fNpoints];
		for (Int_t i = 0; i < fNpoints; i++)
			fX[i] = xaxis->GetBinCenter(i + 1);

		TGraph* pullGr = new TGraph(fNpoints,fX,res);
		pullGr->SetName("dataAdaptiveBinned");
		pad_pull->cd();

		//	pad_pull->SetGridy();
		pullGr->Draw("AP"); //draw markers only
		pullGr->SetTitle("Normalized residuals");
		pullGr->GetXaxis()->SetTitle(xaxis->GetTitle());
		pullGr->GetXaxis()->SetTitleOffset(0.8);
		pullGr->GetXaxis()->SetTitleSize(0.15);
		pullGr->GetXaxis()->SetLabelSize(0.10);
		pullGr->GetXaxis()->SetLimits(
				xaxis->GetBinLowEdge(1),
				xaxis->GetBinLowEdge(fNpoints+1)
		);
		pullGr->GetYaxis()->SetTitle("Deviation [#sigma]");
		pullGr->GetYaxis()->SetTitleOffset(0.4);
		pullGr->GetYaxis()->SetNdivisions(504);
		if(min != max) pullGr->GetYaxis()->SetRangeUser(min,max);
		else pullGr->GetYaxis()->SetRangeUser(-5,5.5);
		pullGr->GetYaxis()->SetTitleSize(0.12);
		pullGr->GetYaxis()->SetLabelSize(0.12);
		pullGr->GetYaxis()->CenterTitle(1);

		TLine* line = new TLine(
				xaxis->GetBinLowEdge(1),
				0.0,
				xaxis->GetBinUpEdge(fNpoints+1),
				0.0
		);
		line->SetLineStyle(2);
		line->Draw();
		gPad->Update();

		gStyle->SetOptTitle(optTitle);
		gStyle->SetOptTitle(optStat);
	}
	return pad;
}

TPad* drawPull(TH1D* hist1, TH1D* hist2, TString drawOption1,
		TString drawOption2, double min=0, double max=0)
{
	std::vector<TH1D*> vHist;
	vHist.push_back(hist1);
	vHist.push_back(hist2);
	std::vector<TString> vOpt;
	vOpt.push_back(drawOption1);
	vOpt.push_back(drawOption2);
	return drawPull(vHist,vOpt,min,max);
}
TPad* drawResidual(std::vector<TH1D*> hist, std::vector<TString> drawOption,
		double min=0, double max=0)
{
	if(hist.size() != drawOption.size() )
		throw std::runtime_error("drawResidual() | Number of histograms and number "
				"of draw options does not match!");

	Int_t optTitle = gStyle->GetOptTitle();
	Int_t optStat = gStyle->GetOptStat();
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);//entries only

	TPad* pad;
	if( !hist.size() ) return 0;
	else if( hist.size() == 1 ) {
		TPad* pad = new TPad();
		hist.at(0)->Draw(drawOption.at(0));
	} else if(hist.size() > 1){
		TPad* pad = new TPad("hist","hist",0.0,0.3,1,1);
		pad->Draw(); pad->SetMargin(0.1,0.05,0.0,0.05);
		TPad* pad_pull = new TPad("pull","pull",0.0,0.0,1,0.3);
		pad_pull->Draw(); pad_pull->SetMargin(0.1,0.05,0.3,0.0);//left-right-bottom-top

		if( hist.at(0)->GetDimension() != 1 || hist.at(1)->GetDimension() != 1 ) {
			std::cout<<"Dimension of histograms larger 1!"<<std::endl;
			return pad;
		}
		if( hist.at(0)->GetNbinsX() != hist.at(1)->GetNbinsX() ) {
			std::cout<<"binning doesnt match"<<std::endl;
			return pad;
		}

		pad->cd();
		hist.at(0)->GetYaxis()->SetRangeUser(0.00000001, hist.at(0)->GetMaximum()*1.1);
		hist.at(0)->GetYaxis()->SetTitleOffset(0.83);

		hist.at(0)->Draw(drawOption.at(0));
		hist.at(1)->Draw(drawOption.at(1));
		for(int i=2; i<hist.size(); ++i)
			hist.at(i)->Draw(drawOption.at(i));

		Double_t chi2;
		Double_t res[hist.at(0)->GetNbinsX()];
		Int_t ndf, igood;
		hist.at(0)->Chi2TestX(hist.at(1),chi2,ndf,igood,"UW",res);

		char chi2Char[60];
		sprintf(chi2Char,"#chi^{2}_{1D}/ndf = %.2f/%d",chi2,ndf);
		TLatex* ltx = new TLatex();
		ltx->DrawLatexNDC(0.2,0.85,chi2Char);//needs to be drawn inside a TPad
		delete ltx;

		TAxis *xaxis = ((TH1*)hist.at(0))->GetXaxis();
		Int_t fNpoints = xaxis->GetNbins();
		//		Double_t fX[fNpoints];
		//		for (Int_t i = 0; i < fNpoints; i++)
		//			fX[i] = xaxis->GetBinCenter(i + 1);

		TH1D* resHist = (TH1D*) hist.at(0)->Clone("resHist");
		resHist->Add(hist.at(1),-1);
		pad_pull->cd();
		resHist->Draw("E"); //draw markers only
		resHist->SetTitle("Residuals");
		resHist->GetXaxis()->SetTitle(xaxis->GetTitle());
		resHist->GetXaxis()->SetTitleOffset(0.8);
		resHist->GetXaxis()->SetTitleSize(0.15);
		resHist->GetXaxis()->SetLabelSize(0.10);
		resHist->GetXaxis()->SetLimits(
				xaxis->GetBinLowEdge(1),
				xaxis->GetBinLowEdge(fNpoints+1)
		);
		resHist->GetYaxis()->SetTitle("Deviation");
		resHist->GetYaxis()->SetTitleOffset(0.4);
		resHist->GetYaxis()->SetNdivisions(504);
		if(min != max) resHist->GetYaxis()->SetRangeUser(min,max);
		resHist->GetYaxis()->SetTitleSize(0.12);
		resHist->GetYaxis()->SetLabelSize(0.12);
		resHist->GetYaxis()->CenterTitle(1);

		TLine* line = new TLine(
				xaxis->GetBinLowEdge(1),
				0.0,
				xaxis->GetBinUpEdge(fNpoints+1),
				0.0
		);
		line->SetLineStyle(2);
		line->Draw();
		gPad->Update();

		gStyle->SetOptTitle(optTitle);
		gStyle->SetOptTitle(optStat);
	}
	return pad;
}

TPad* drawResidual(TH1D* hist1, TH1D* hist2, TString drawOption1,
		TString drawOption2, double min=0, double max=0)
{
	std::vector<TH1D*> vHist;
	vHist.push_back(hist1);
	vHist.push_back(hist2);
	std::vector<TString> vOpt;
	vOpt.push_back(drawOption1);
	vOpt.push_back(drawOption2);
	return drawResidual(vHist,vOpt,min,max);
}

void getTH2PolyChi2(TH2Poly* hist1, TH2Poly* hist2, double& chi2, int& ndf, int& igood)
{
	if(hist1->GetBins()->GetEntries() != hist2->GetBins()->GetEntries() ){
		std::cout<<"binning doesnt match"<<std::endl;
		return;
	}
	//	std::cout<<hist1->GetEntries()<<std::endl;
	//	std::cout<<hist1->GetIntegral()<<std::endl;
	//	std::cout<<hist1->GetSumOfWeights()<<std::endl;
	chi2 = 0;
	ndf = 0;
	for(int bin=1; bin<=hist1->GetBins()->GetEntries(); bin++) {
		double c1 = hist1->GetBinContent(bin);
		double c2 = hist2->GetBinContent(bin);
		double c1Err = hist1->GetBinError(bin);
		double c2Err = hist2->GetBinError(bin);
		if(c1>0 || c2>0) {
			chi2+=(c1-c2)*(c1-c2)/(c1Err*c1Err+c2Err*c2Err);
			//			chi2+=(c1-c2)*(c1-c2)/(c1Err*c1Err);
			//			std::cout<<chi2<<std::endl;
			ndf++;
		}
	}
	//	std::cout<<"Calculating TH2Poly chi2="<<chi2<<" ndf="<<ndf<<std::endl;
	return;
}

TH2Poly* getTH2PolyPull(TH2Poly* hist1, TH2Poly* hist2, TString name){
	if(hist1->GetBins()->GetEntries() != hist2->GetBins()->GetEntries() )
	{
		std::cout<<"binning doesnt match"<<std::endl;
		return 0;
	}
	TH2Poly* resHist = (TH2Poly*)hist1->Clone(name+hist1->GetName());
	resHist->Reset("");

	hist2->Scale(hist1->Integral()/hist2->Integral());
	double limit=0;
	double integral=0;
	int nBins=0;
	for(int bin=1; bin<=hist1->GetBins()->GetEntries(); bin++) {
		double c1 = hist1->GetBinContent(bin);
		double c2 = hist2->GetBinContent(bin);
		double c1Err = hist1->GetBinError(bin);
		double c2Err = hist2->GetBinError(bin);
		if(c1>0 && c2>0) {
			double pull = (c1-c2)/sqrt(c1Err*c1Err+c2Err*c2Err);
			integral+=pull;
			nBins++;
			resHist->SetBinContent(bin,pull);
			if(std::fabs(pull)>limit) limit=std::abs(pull);
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
	std::cout<<"Creating pull histogram for TH2Poly! integral="<<integral<<" nBins="<<nBins<<std::endl;

	return resHist;
}




TH2Poly* adaptiveBinning(UInt_t dataSize, UInt_t dataDim, Double_t* data, UInt_t nBins = 100)
{
	UInt_t size = UInt_t(dataSize/nBins)*nBins;//size should be a multiple of nBins
	if( size==0 ) { //return empty histogram
		TH2Poly *h2p = new TH2Poly();
		h2p->AddBin(0, 1, 0, 1);
		return h2p;
	}

	TKDTreeBinning kdBins(size, dataDim, data, nBins);

	UInt_t nbins = kdBins.GetNBins();
	UInt_t dim   = kdBins.GetDim();

	const Double_t* binsMinEdges = kdBins.GetBinsMinEdges();
	const Double_t* binsMaxEdges = kdBins.GetBinsMaxEdges();

	TH2Poly* h2pol = new TH2Poly("h2PolyBinTest", "",
			kdBins.GetDataMin(0), kdBins.GetDataMax(0), kdBins.GetDataMin(1), kdBins.GetDataMax(1));

	for (UInt_t i = 0; i < nbins; ++i) {
		UInt_t edgeDim = i * dim;
		h2pol->AddBin(binsMinEdges[edgeDim], binsMinEdges[edgeDim + 1],
				binsMaxEdges[edgeDim], binsMaxEdges[edgeDim + 1]);
	}

	return h2pol;
}
#endif
