
void xy_graph(vector<double> x, vector <double> y);

void benchmarkMacro (){
    int n_file = 20;
    TChain* chain = new TChain("SimulatedTree");
    for (auto i=0; i<n_file; i++){
        string filename = "root/simulation/benchmark" + to_string(i+1) + ".root";
        chain->Add(filename.c_str());
    }

    TSamuraiIdealData* FDC1 = new TSamuraiIdealData();
    TSamuraiIdealData* FDC2 = new TSamuraiIdealData();
    TSamuraiIdealData* Magnet = new TSamuraiIdealData();

    chain->SetBranchAddress("IdealDataFDC1", &FDC1);
    chain->SetBranchAddress("IdealDataFDC2", &FDC2);
    chain->SetBranchAddress("IdealDataMagnet", &Magnet);

    vector <double> X, Y;

    int entries = chain->GetEntries();
    for (int i=0; i<entries; i++){
        chain->GetEntry(i);
        for (auto j=0;j<Magnet->GetMult(); j++){
            X.push_back(-Magnet->GetPosX(j));
            //Y.push_back(FDC2->GetPosY(j));
            Y.push_back(Magnet->GetPosZ(j));
            //Y.push_back(FDC2->GetBrho(j));
        }

    }

    double ylim = 20, xa=-3812, xb=-3794;
    TH2D* h2d = new TH2D("h2d", "X:Y", 100, xa, xb, 100, -ylim, ylim);
    //TH2D* h2d = new TH2D("h2d", "X:Brho", 100, -3812, -3794, 100, 6.06, 6.12);
    
    for (int i=0; i<X.size(); i++){
        h2d->Fill(X[i], Y[i]);
    }
    h2d->Draw("colz");
    xy_graph(X, Y);

}

//-------------------------------------------------------------

void t_hist(vector<double> t){

    TH1D* hist = new TH1D("h1", "h1 title", 50, 0, 1);
    for (unsigned int i = 0; i < t.size(); i++) hist->Fill(t[i]);

    TF1* tf = new TF1("tf", "[0]", 0, 1);
    tf->SetParameters(1, 1);
    tf->SetParNames("A", "B");
    hist->Fit(tf);
    
    cout << "chi square : " << tf->GetChisquare() << endl;
    cout << "prob : " << tf->GetProb() << endl;

    TCanvas* tc = new TCanvas();
    hist->GetYaxis()->SetRangeUser(0, 1.5*t.size()/50);
    hist->Draw();
    tc->Print("hist.pdf");

    return;
}

void xy_graph(vector <double> x, vector <double> y){

    TGraph* tg = new TGraph(x.size(), &x[0], &y[0]);

    //TF1* tf = new TF1("tf", "pol2", 0, 1);
    //tf->SetParameters(1, 1, 1);
    //tf->SetParNames("A", "B", "C");
    //tg->Fit(tf);

    tg->SetMarkerStyle('+'); //draws x

    //TCanvas* tc = new TCanvas();
    tg->Draw("AP");
    //tc->Print("graph.pdf");

    return;
}

void tx_graperrors(vector <double> t, vector <double> x, vector <double> et, vector <double> ex){

    TGraphErrors* tg = new TGraphErrors(t.size(), &t[0], &x[0], &et[0], &ex[0]);

    TF1* tf = new TF1("tf", "pol2", 0, 1);
    tf->SetParameters(1, 1, 1);
    tf->SetParNames("A", "B", "C");
    tg->Fit(tf);

    cout << "chi square : " << tf->GetChisquare() << endl;
    cout << "prob : " << tf->GetProb() << endl;

    TCanvas* tc = new TCanvas();
    tg->Draw("AP");
    tc->Print("grapherrors.pdf");
}
