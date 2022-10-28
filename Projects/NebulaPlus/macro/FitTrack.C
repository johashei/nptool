TChain* chain=NULL ;

////////////////////////////////////////////////////////////////////////////////
void LoadChain(){
    chain = new TChain("PhysicsTree");
    //chain->Add("./root/analysis/testmult16.root");
    chain->Add("./root/analysis/run_38_0001.root");
}
double fpol1(double* x, double* par){
    double result = par[0]+par[1]*x[0];
    return result;
}

////////////////////////////////////////////////////////////////////////////////
void FitTrack(){
    LoadChain();

    TNebulaPlusPhysics* m_NebulaPlus = new TNebulaPlusPhysics();
    chain->SetBranchAddress("NebulaPlus", &m_NebulaPlus);
    Long64_t nentries = chain->GetEntries();
    TH1F* h = new TH1F("h","h",300,0,300);
    TH1F* hTheta = new TH1F("hTheta","hTheta",100,0,100);
    TH1F** hdT = new TH1F*[90];
    for ( int i=0; i<90; i++) {
        hdT[i]=new TH1F(Form("hdT%d",i+1),Form("TOF%d - TOF%d:TOF(ns)",i+1,i+2),400,-10,10);
    }

    TCanvas* c = new TCanvas("c","c",0,0,700,700);
    c->cd();
    gStyle->SetOptFit(111111111);
    h->GetXaxis()->SetRangeUser(-300,300);
    h->GetYaxis()->SetRangeUser(-300,300);
    double Theta;

    for(Long64_t i = 0; i < nentries; i++){
        chain->GetEntry(i);
        int mult= m_NebulaPlus->DetectorNumber.size();
        
        /// build flag to select single layer event ///
        double tmp=-99;
        bool multilayer=false;
        for(int j=0;j<mult;j++){
            if(j==0) tmp = m_NebulaPlus->PosZ[j];
            if(m_NebulaPlus->PosZ[j]!=tmp) multilayer=true; 
        }

        if(!multilayer){
            if(mult>6){

                ///////////////////////////////////////////////
                //// Track construction and Fit in Y vs X /////
                ///////////////////////////////////////////////
                TGraph* g = new TGraph(
                        mult,
                        &(m_NebulaPlus->PosX[0]),
                        &(m_NebulaPlus->PosY[0])
                        );
                g->SetMarkerStyle(3);
                h->Draw("");
                g->Draw("Psame");

                TF1* toto = new TF1("f1",fpol1,0,300,2);
                g->Fit("f1","QR");
                double offset = toto->GetParameter(0);
                double slope = toto->GetParameter(1);

                //std::cout <<  offset << std::endl;;
                //std::cout <<  slope << std::endl;;

                TVector2 a(0,offset);
                TVector2 b(1,offset+slope);
                TVector2 c = b-a;
                Theta = fabs(TMath::ATan(c.Y()/c.X())*180./TMath::Pi());
                std::cout << Theta << std::endl;
                hTheta->Fill(Theta);
                gPad->WaitPrimitive();
                
                ///////////////////////////////////////////////
                //// Time analysis between adjacent bars //////
                ///////////////////////////////////////////////
                
                vector<int> ID(mult);     
                vector<int> IDsort(mult);     
                vector<int> idx(mult);     
                vector<double> T(mult);     
                double dT=0;     

                ID = m_NebulaPlus->DetectorNumber;

                //std::cout << "---------" << std::endl;
                //std::cout << ID[0] << std::endl;
                //std::cout << m_NebulaPlus->DetectorNumber[0] << std::endl;
                TMath::Sort(mult,&ID[0],&idx[0],false);

                //std::cout << "---------" << std::endl;

                for(int k=0;k<mult;k++){
                    T[k] =  m_NebulaPlus->TOF[idx[k]];
                    IDsort[k] =  m_NebulaPlus->DetectorNumber[idx[k]];
                    //if(k=0) refT = T[0]; refID = ID[0];
                    if(k>0 && IDsort[k-1]==IDsort[k]-1 ){
                        // if(T[k]>T[k-1]) dT = T[k]-T[k-1];    
                        // if(T[k]<T[k-1]) dT = T[k-1]-T[k];    

                        dT = T[k]-T[k-1];    
                        hdT[IDsort[k-1]]->Fill(dT);

                        //std::cout <<IDsort[k]<< "\t" << hdT[IDsort[k-1]]->GetName() << "\t" << dT<< std::endl;
                    }
                }



            }
        }

    }

    c->cd();
    hTheta->Draw();
    c->Print("TOF.pdf[");
    TGraph *res1 = new TGraph(90);
    TGraph *res2 = new TGraph(90);
    res1->SetName("res1");
    res2->SetName("res2");
    for (int i=0; i<90; i++){
        TF1* f_gauss = new TF1("f_gauss","gaus(0)+gaus(3)",-10,10);
        double par[6] = {hdT[i]->GetMaximum(), hdT[i]->GetMean()-.5, 0.5,hdT[i]->GetMaximum(), hdT[i]->GetMean()+0.5 , 0.5 };
        f_gauss->SetParameters(par);
        f_gauss->SetNpx(1000);
        hdT[i]->Fit(f_gauss,"Q");
        c->Print("TOF.pdf");
        res1->SetPoint(i,i+1,TMath::Abs(f_gauss->GetParameter(2)));
        res2->SetPoint(i,i+1,TMath::Abs(f_gauss->GetParameter(5)));

    }
    res1->SetMarkerStyle(kStar);
    res1->SetMarkerColor(kRed);


    res2->SetMarkerStyle(kStar);
    res2->SetMarkerColor(kBlue);
    TMultiGraph* mg = new TMultiGraph();
    mg->Add(res1,"p");
    mg->Add(res2,"p");
    mg->SetTitle(" ; ID ; sigma value");
    mg->Draw("a");


    c->Print("TOF.pdf");
    c->Print("TOF.pdf]");





}
