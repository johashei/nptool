
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

    //vector<string> paths ={"spm1000/1cm_E0/", "spm1000/1cm_E480/", "spm1000/5cm_E480/", 
    //                     "spm10000/1cm_E0/", "spm10000/1cm_E480/", "spm10000/5cm_E480/"};
    //vector<string> paths ={"spm1000/1cm_E0/", "spm1000/1cm_E480/", "spm1000/5cm_E480/"};
    vector<string> paths ={"spm10000/1cm_E0/", "spm10000/1cm_E480/", "spm10000/5cm_E480/"};
    for(auto p:paths) Brhos(p);

    //Brhos("spm1000/1cm_E0/");

}

void Brhos(string path_to_file){

    cout << "Read from " << path_to_file << endl;

    string filenameIn = path_to_file + "testanalysis.root";
    string filenameOut = path_to_file + "brhos.root";
    string fieldmap_file = "../../field_map/3T.table.bin";    

    TFile* fileIn = new TFile(filenameIn.c_str(), "read");
    TTree* input = (TTree*)fileIn->Get("SimulatedTree");
    TFile* fileOut = new TFile(filenameOut.c_str(), "recreate");
    TTree* output = new TTree("BrhoTree", "brhos");

    SamuraiFieldMap* Map = new SamuraiFieldMap();
    double mag_angle = 30*deg, fdc2_angle=-(59.930-30)*deg;
    Map->LoadMap(mag_angle, fieldmap_file, 10);
    Map->SetFDC2R((5000-438)*mm); //z pos - half thickness of FDC2
    Map->SetFDC2Angle(fdc2_angle);
    Map->SetStepLimit(1e6);

    TSamuraiIdealData* FDC1 = new TSamuraiIdealData();
    TSamuraiIdealData* FDC2 = new TSamuraiIdealData();
    TSamuraiIdealData* Magnet = new TSamuraiIdealData();

    input->SetBranchAddress("IdealDataFDC1", &FDC1);
    input->SetBranchAddress("IdealDataFDC2", &FDC2);
    input->SetBranchAddress("IdealDataMagnet", &Magnet);


    const vector <int> time = {
        10, 20, 40, 70, 100, 200, 300, 400, 600, 800, 
        1000, 2000, 3000, 4000, 6000, 8000, 10000, 20000, 
        30000, 40000, 60000, 80000, 100000, 150000, 200000};
    //const vector <int> time = {100};
    double br_recs [time.size()];

    for (auto i=0; i<time.size(); i++){
        output->Branch( ("Brho"+to_string(time[i])).c_str(), &br_recs[i], "brho/D");
    }

    int entries = input->GetEntries();

    TVector3 Mag_Pos (0,0,10000*mm);
    TVector3 Pos_FDC1, Dir_FDC1, Pos_FDC2, Dir_FDC2;
    for (int i=0; i<entries; i++){
        cout << path_to_file << i << " ";
        for (auto j=0; j<time.size(); j++) br_recs[j]=-1;
        input->GetEntry(i);
        if (FDC2->GetMult() == 1 && FDC1->GetMult() == 1){
            Pos_FDC1 = TVector3(FDC1->GetPosX(0), FDC1->GetPosY(0), FDC1->GetPosZ(0)) - Mag_Pos;
            Pos_FDC2 = TVector3(FDC2->GetPosX(0), FDC2->GetPosY(0), FDC2->GetPosZ(0)) - Mag_Pos;
            Dir_FDC1.SetMagThetaPhi(FDC1->GetMomMag(0), FDC1->GetMomTheta(0), FDC1->GetMomPhi(0));
            Dir_FDC2.SetMagThetaPhi(FDC2->GetMomMag(0), FDC2->GetMomTheta(0), FDC2->GetMomPhi(0));
            Dir_FDC1.Unit();
            Dir_FDC2.Unit();
            Dir_FDC2.RotateY(mag_angle-fdc2_angle);
            Pos_FDC2.RotateY(mag_angle-fdc2_angle);
            
            for (auto j=0; j<time.size(); j++){
                Map->SetTimeIntervalSize(time[j]*picosecond);
                br_recs[j] = Map->FindBrho(Pos_FDC1, Dir_FDC1, Pos_FDC2, Dir_FDC2);
            }
        }
        output->Fill();
    }
    fileIn->Close();
    fileOut->Write();
    fileOut->Close();
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