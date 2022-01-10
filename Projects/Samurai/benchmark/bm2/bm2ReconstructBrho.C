void Brhos(string path_to_file);

void bm2ReconstructBrho (){

    //vector<string> paths ={"field_10mm/spm1000/1cm_E0/", "field_10mm/spm1000/1cm_E480/", "field_10mm/spm1000/5cm_E480/", 
    //                     "field_10mm/spm10000/1cm_E0/", "field_10mm/spm10000/1cm_E480/", "field_10mm/spm10000/5cm_E480/"};
    //vector<string> paths ={"field_10mm/1x5y/", "field_20mm/"};
    //for(auto p:paths) Brhos(p);

    Brhos("resolution/1cm_E480/");
    //Brhos("field_10mm/5x1y/");
    //Brhos("field_20mm/");
    //Brhos("field_10mm/spm1000/1cm_E0/");

    //vector<string> paths ={"field_10mm/1x5y/", "field_20mm/"};
    //for(auto p:paths) Brhos(p);

}

void Brhos(string path_to_file){

    cout << "Read from " << path_to_file << endl;

    string filenameIn = path_to_file + "testanalysis.root";
    string filenameOut = path_to_file + "brhos.root";
    //string fieldmap_file = "../../field_map/3T_20mm.table.bin";
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
    double fx1, fx2, fy1, fy2;
    static const double Res_FDC1 = 110*um, Res_FDC2 = 120*um;
    static TRandom2* Rand_Generator = new TRandom2 ();
    for (int i=0; i<entries; i++){
        cout << i << " " << path_to_file  << " ";
        for (auto j=0; j<time.size(); j++) br_recs[j]=-1;
        input->GetEntry(i);
        if (FDC2->GetMult() == 1 && FDC1->GetMult() == 1){
            fx1 = Rand_Generator->Gaus(FDC1->GetPosX(0), Res_FDC1);
            fy1 = Rand_Generator->Gaus(FDC1->GetPosY(0), Res_FDC1);
            fx2 = Rand_Generator->Gaus(FDC2->GetPosX(0), Res_FDC2);
            fy2 = Rand_Generator->Gaus(FDC2->GetPosY(0), Res_FDC2);
            Pos_FDC1 = TVector3(fx1, fy1, FDC1->GetPosZ(0)) - Mag_Pos;
            Pos_FDC2 = TVector3(fx2, fy2, FDC2->GetPosZ(0)) - Mag_Pos;
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
