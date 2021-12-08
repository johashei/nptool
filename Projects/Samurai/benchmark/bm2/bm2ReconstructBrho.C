
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


void Brhos(string path_to_file);
void findbrho(string path_to_file, int time_int, SamuraiFieldMap* Map, double mag_angle, double fdc2_angle,
            const vector<double>& Brho, const vector<TVector3>& Pos_FDC1,
            const vector<TVector3>& Pos_FDC2, const vector<TVector3>& Dir_FDC1,
            const vector<TVector3>& Dir_FDC2);


void bm2ReconstructBrho (){

   vector<string> paths ={"spm1000/1cm_E0/", "spm1000/1cm_E480/", "spm1000/5cm_E480/", 
                        "spm10000/1cm_E0/", "spm10000/1cm_E480/", "spm10000/5cm_E480/"};
   for(auto p:paths) Brhos(p);
    
    /*
    for (auto i =0; i<Br_rec.size(); i++) if(Br_rec[i] < 6) cout << i << " " << Br_rec[i] << endl;
    TH1D* hist1d = new TH1D("brho", "Brho", 200, 6, 7);
    for (int i = 0; i < Br_rec.size(); i++) hist1d->Fill(Br_rec[i]);
    TCanvas* tc = new TCanvas();
    tc->SetLogy();
    hist1d->Draw();*/
/*
    vector<TVector3> track1 = Map->Propagate(6.26709, Pos_FDC1[7839], Dir_FDC1[7839], true);
    vector<TVector3> track2 = Map->Propagate(Br_rec[839], Pos_FDC1[7839], Dir_FDC1[7839], true);
    vector<TVector3> track3 = Map->Propagate(1, Pos_FDC1[7839], Dir_FDC1[7839], true);

    vector<double> x1, z1, x2, z2, x3, z3;
    for (int i=0; i<track1.size(); i++){
        x1.push_back(-track1[i].X());
        z1.push_back(track1[i].Z());
    }
    for (int i=0; i<track2.size(); i++){
        x2.push_back(-track2[i].X());
        z2.push_back(track2[i].Z());
    }
    for (int i=0; i<track3.size(); i++){
        x3.push_back(-track3[i].X());
        z3.push_back(track3[i].Z());
    }

    TGraph* tg1 = new TGraph(x1.size(), &x1[0], &z1[0]);
    TGraph* tg2 = new TGraph(x2.size(), &x2[0], &z2[0]);
    TGraph* tg3 = new TGraph(x3.size(), &x3[0], &z3[0]);
    TMarker* mk = new TMarker(-Pos_FDC2[7839].X(), Pos_FDC2[7839].Z(), 29);
    tg1->SetMarkerStyle(2);
    tg2->SetMarkerStyle(3);
    tg3->SetMarkerStyle(4);

    tg1->Draw("APL");
    tg2->Draw("PLsame");
    tg3->Draw("PLsame");
    mk->Draw();
    */
    
   // cout << "Real position FDC1: " << "\t" << Pos_FDC1[0].X() << " " << Pos_FDC1[0].Y() << " " << Pos_FDC1[0].Z() << endl;
   // cout << "Real direction FDC1: " << "\t" << Dir_FDC1[0].X() << " " << Dir_FDC1[0].Y() << " " << Dir_FDC1[0].Z() << endl;
   // cout << "Real position FDC2: " << "\t" << Pos_FDC2[0].X() << " " << Pos_FDC2[0].Y() << " " << Pos_FDC2[0].Z() << endl;
   // cout << "Real direction FDC2: " << "\t" << Dir_FDC2[0].X() << " " << Dir_FDC2[0].Y() << " " << Dir_FDC2[0].Z() << endl;
   // cout << "I found Brho = " << br << endl;
   // cout << "Real Brho = " << Brho[0] << endl;

    

    //TH1D* h1d = new TH1D("h1d", "Brho at FDC2", 100, 5, 8);
    //for (auto i=0; i<Brho.size(); i++) h1d->Fill(Brho[i]);
    //h1d->Draw();
    //tc->Print((path+set+"histoBrho_FDC2.pdf").c_str());

   // cout << "Brho\n";
   // cout << "mean = " << mean(Brho) << endl;
   // cout << "std dev = " << sqrt(var(Brho)) << endl;
   // 
   // cout << "Brho Reconstructed\n";
   // cout << "mean = " << mean(Br_rec) << endl;
   // cout << "std dev = " << sqrt(var(Br_rec)) << endl;

}

void Brhos(string path_to_file){
    /*string path = "";
    string set = "spm1000/1cm_E0/";
    string filename = path + set + "testanalysis.root";*/

    //cout << "Read from " << path + set << endl;
    cout <<"Read from "<<path_to_file << endl;

    string filename = path_to_file + "testanalysis.root";
    

    string fieldmap_file = "../../field_map/3T.table.bin";    
    //string fieldmap_file = path_to_fieldmap + "3T.table.bin";
    TFile* file = new TFile(filename.c_str(), "read");
    TTree* tree = (TTree*)file->Get("SimulatedTree");

    SamuraiFieldMap* Map = new SamuraiFieldMap();
    TSamuraiIdealData* FDC1 = new TSamuraiIdealData();
    TSamuraiIdealData* FDC2 = new TSamuraiIdealData();
    TSamuraiIdealData* Magnet = new TSamuraiIdealData();

    tree->SetBranchAddress("IdealDataFDC1", &FDC1);
    tree->SetBranchAddress("IdealDataFDC2", &FDC2);
    tree->SetBranchAddress("IdealDataMagnet", &Magnet);

    vector <double> Brho;
    vector <TVector3> Pos_FDC1, Dir_FDC1, Pos_FDC2, Dir_FDC2;

    int entries = tree->GetEntries();

    TVector3 Mag_Pos (0,0,10000*mm);
    TVector3 dir1, dir2;
    for (int i=0; i<entries; i++){
        tree->GetEntry(i);
        if (FDC1->GetMult() == FDC2->GetMult()){
            for (auto j=0;j<FDC2->GetMult(); j++){
                Brho.push_back(FDC2->GetBrho(j));
                Pos_FDC1.push_back(TVector3(FDC1->GetPosX(j), FDC1->GetPosY(j), FDC1->GetPosZ(j) ) - Mag_Pos);
                Pos_FDC2.push_back(TVector3(FDC2->GetPosX(j), FDC2->GetPosY(j), FDC2->GetPosZ(j) ) - Mag_Pos);
                dir1.SetMagThetaPhi(FDC1->GetMomMag(j), FDC1->GetMomTheta(j), FDC1->GetMomPhi(j));
                Dir_FDC1.push_back(dir1.Unit());
                dir2.SetMagThetaPhi(FDC2->GetMomMag(i), FDC2->GetMomTheta(j), FDC2->GetMomPhi(j));
                Dir_FDC2.push_back(dir2.Unit());
            }
        }
    }
    //int time_int = 1000;
    double mag_angle = 30*deg, fdc2_angle=-(59.930-30)*deg;
    Map->LoadMap(mag_angle, fieldmap_file, 10);
    Map->SetFDC2R((5000-438)*mm); //z pos - half thickness of FDC2
    Map->SetFDC2Angle(fdc2_angle);
    Map->SetStepLimit(1e5);

    vector <int> time {10000,3000,1000,300,100};
    //vector <int> time {10000,3000};
    for(auto t:time) findbrho(path_to_file, t, Map, mag_angle, 
                    fdc2_angle, Brho, Pos_FDC1, Pos_FDC2, Dir_FDC1, Dir_FDC2);
////////////////////////////////////////////////////////////
    

}

void findbrho(string path_to_file, int time_int, SamuraiFieldMap* Map, double mag_angle, double fdc2_angle,
            const vector<double>& Brho, const vector<TVector3>& Pos_FDC1,
            const vector<TVector3>& Pos_FDC2, const vector<TVector3>& Dir_FDC1,
            const vector<TVector3>& Dir_FDC2)
{
    Map->SetTimeIntervalSize(time_int*picosecond);

    TVector3 posF2, dirF2;
    vector <double> Br_rec;
    for (auto i=0; i<Brho.size();i++){
        cout<<i<<" ";
        posF2 = Pos_FDC2[i];
        posF2.RotateY(mag_angle-fdc2_angle);
        dirF2 = Dir_FDC2[i];
        dirF2.RotateY(mag_angle-fdc2_angle);
        double brs = Map->FindBrho(Pos_FDC1[i], Dir_FDC1[i], posF2, dirF2);
        Br_rec.push_back(brs);
    }
    double brho;
    TFile* file2 = new TFile((path_to_file + "brho_" + to_string(time_int) + "ps.root").c_str(), "recreate");
    TTree* brhos = new TTree("tree", "brhos");
    brhos->Branch("Brho", &brho, "brho/D");
    for(auto i=0; i<Br_rec.size(); i++) {
        brho = Br_rec[i];
        brhos->Fill();
    }
    file2->Write();
    file2->Close();
}