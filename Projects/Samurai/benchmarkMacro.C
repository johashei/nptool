
void xy_graph(vector<double> x, vector <double> y);

double mean(vector<double>V){
    double s = 0;
    for (auto x:V) s += x;
    return s / V.size();
}
double var(vector<double>V){
    double s = 0, m=mean(V);
    for (auto x:V) s += (x-m)*(x-m);
    return s / (V.size()-1);
}


void benchmarkMacro (){
    int n_file = 20;
    string path = "root/simulation/bm1/";
    string set = "EO/spm3/";
    cout << "Read from " << path + set << endl;
    TChain* chain = new TChain("SimulatedTree");
    for (auto i=0; i<n_file; i++){
        string filename = path + set + "benchmark" + to_string(i+1) + ".root";
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
        for (auto j=0;j<FDC2->GetMult(); j++){
            X.push_back(FDC2->GetPosX(j));
            Y.push_back(FDC2->GetPosY(j));
            //Y.push_back(FDC2->GetPosZ(j));
            //Y.push_back(FDC2->GetBrho(j));
        }

    }

    TCanvas* tc = new TCanvas();

    TH1D* h1d = new TH1D("h1d", "X at FDC2", 100, -3812, -3794);
    for (auto i=0; i<X.size(); i++) h1d->Fill(X[i]);
    h1d->Draw();
    tc->Print((path+set+"histoX_FDC2.pdf").c_str());

    TH2D* h2d = new TH2D("h2d", "X:Y", 100, -3812, -3794, 100, -20, 20);
    for (int i=0; i<X.size(); i++) h2d->Fill(X[i], Y[i]);
    h2d->Draw("colz");
    tc->Print((path+set+"histoXY_FDC2.pdf").c_str());

    cout << "X\n";
    cout << "mean = " << mean(X) << endl;
    cout << "std dev = " << sqrt(var(X)) << endl;
    cout << "Y\n";
    cout << "mean = " << mean(Y) << endl;
    cout << "std dev = " << sqrt(var(Y)) << endl;
    

}

//-------------------------------------------------------------


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
