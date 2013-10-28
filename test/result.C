{
  gROOT->Reset();
  gStyle->SetOptTitle(0);
  TCanvas *cc = new TCanvas("cc","cc",0,0,1024,768);
  cc->Divide(1,2);

  std::vector<long> iteration;

  std::vector<double> evaluation0;

  std::vector<double> evaluation1;

  iteration.push_back(0);
  evaluation0.push_back(3.25851e+07);
  evaluation1.push_back(3.35735e+07);

  iteration.push_back(1000);
  evaluation0.push_back(32294.3);
  evaluation1.push_back(32296.6);

  iteration.push_back(2000);
  evaluation0.push_back(30917.4);
  evaluation1.push_back(30917.4);

  iteration.push_back(3000);
  evaluation0.push_back(26124);
  evaluation1.push_back(26124);

  iteration.push_back(4000);
  evaluation0.push_back(25877.8);
  evaluation1.push_back(25877.8);

  iteration.push_back(5000);
  evaluation0.push_back(25503.8);
  evaluation1.push_back(25503.8);

  iteration.push_back(6000);
  evaluation0.push_back(25503.8);
  evaluation1.push_back(25503.8);

  iteration.push_back(7000);
  evaluation0.push_back(25063.4);
  evaluation1.push_back(25063.4);

  iteration.push_back(8000);
  evaluation0.push_back(24336.9);
  evaluation1.push_back(24336.9);

  iteration.push_back(9000);
  evaluation0.push_back(23544.3);
  evaluation1.push_back(23544.3);

  iteration.push_back(10000);
  evaluation0.push_back(19749.1);
  evaluation1.push_back(19762.9);

  // Transfer the vectors into arrays
  double iteration_arr[iteration.size()];
  double evaluation0_arr[evaluation0.size()];

  for(std::size_t i=0; i<iteration.size(); i++) {
     iteration_arr[i] = (double)iteration[i];     evaluation0_arr[i] = evaluation0[i];
  }

  // Create a TGraph object
  TGraph *evGraph0 = new TGraph(evaluation0.size(), iteration_arr, evaluation0_arr);
  // Set the axis titles
  evGraph0->GetXaxis()->SetTitle("Iteration");
  evGraph0->GetYaxis()->SetTitleOffset(1.1);
  evGraph0->GetYaxis()->SetTitle("Fitness");

  double evaluation1_arr[evaluation1.size()];

  for(std::size_t i=0; i<iteration.size(); i++) {
     iteration_arr[i] = (double)iteration[i];     evaluation1_arr[i] = evaluation1[i];
  }

  // Create a TGraph object
  TGraph *evGraph1 = new TGraph(evaluation1.size(), iteration_arr, evaluation1_arr);
  // Set the axis titles
  evGraph1->GetXaxis()->SetTitle("Iteration");
  evGraph1->GetYaxis()->SetTitleOffset(1.1);
  evGraph1->GetYaxis()->SetTitle("Fitness");

  // Do the actual drawing
  cc->cd(1);
  evGraph0->Draw("AP");
  cc->cd(2);
  evGraph1->Draw("AP");
  cc->cd();
}
