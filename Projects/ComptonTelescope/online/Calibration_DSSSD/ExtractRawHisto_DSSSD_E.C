////////////////////////////////////////////////////////////////////////////
// This macro takes a converted ROOT file and                             //
// create an histogram for each strip (X and Y) of the DSSSD array filled //
// with the energy. The histograms are dumped in an output ROOT file.     //
//                                                                        //
// This is especially usefull for calibration purposes when there is no   //
// need to work directly on the TChain                                    //
////////////////////////////////////////////////////////////////////////////

// C++ headers
#include <iostream>
#include <fstream>
using namespace std;

// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH2.h"

// custom headers
#include "../src/DecodeD.h"
#include "../src/DecodeD.cpp"

#define NBDETECTORS	1	
#define	NBSTRIPS	32

void ExtractRawHisto_DSSSD_E(const char* filename = "20200128_10h44_bi207_conv")
{

   // open the output ROOT file
//   TString outFileName = "./Histograms/";
   TString outFileName = filename;
   outFileName += "_RawDSSSDHistos.root";
   cout << "output file: " << outFileName << endl;
   TFile *outFile = new TFile(outFileName, "recreate");

   // prepare output histograms
   // individual strips
   TH1F* hFrontEnergy[NBDETECTORS][NBSTRIPS];
   TH1F* hBackEnergy[NBDETECTORS][NBSTRIPS];
   // individual detector
   TH2D* hFront[NBDETECTORS];
   TH2D* hBack[NBDETECTORS];
   for (Int_t i = 0; i < NBDETECTORS; i++) {
     // detector Front
     TString hnameF  = Form("h_D%i_FRONT_ENERGY", i+1);
     TString htitleF = Form("D%i_FRONT_ENERGY", i+1);
     hFront[i]       = new TH2D(hnameF, htitleF, NBSTRIPS, 0, NBSTRIPS, 1024, 0, 1024);
     // detector Back
     TString hnameB  = Form("h_D%i_BACK_ENERGY", i+1);
     TString htitleB = Form("D%i_BACK_ENERGY", i+1);
     hBack[i]        = new TH2D(hnameB, htitleB, NBSTRIPS, 0, NBSTRIPS, 1024, 0, 1024);
     for (Int_t j = 0; j < NBSTRIPS; j++) {
       // strips Front
       TString hnameFE     = Form("h_D%d_FRONT_E%d", i+1, j+1);
       TString htitleFE    = Form("D%d_FRONT_E%d", i+1, j+1);
       hFrontEnergy[i][j] = new TH1F(hnameFE, htitleFE, 1024, 0, 1024);
       // strips Back
       TString hnameBE     = Form("h_D%d_BACK_E%d", i+1, j+1);
       TString htitleBE    = Form("D%d_BACK_E%d", i+1, j+1);
       hBackEnergy[i][j] = new TH1F(hnameBE, htitleBE, 1024, 0, 1024);
     }
   }

   // Instantiates DecodeD object reading DSSSD(s) data flux
   cout << "Reading data" << endl;
   DecodeD* DD = new DecodeD(true);
   newframe_t* event;

	 // run to analyse
	 DD -> setTree(Form("../../data/%s.root",filename));
   int dlength = DD -> getLength();

   while (DD -> getCursor() < dlength) {
     DD -> decodeEvent();
     event = DD -> getEvent();
     if (DD -> getCursor() % 10000 == 0) {
       cout << "\rEntry " << DD -> getCursor();
       cout << flush;
     }

     
     for (int i = 0; i < 2; i++) { // 2 faces
       for (int j = 0; j < NBDETECTORS; j++) { // detectors
         if (event->chip_data[i][j]) { // Test if data is present
           switch (i) {
             case 0: //Assuming 0 is front - to be checked
               for (int k = 0; k < NBSTRIPS; k++) { // strips         
                 hFrontEnergy[j][k]->Fill(event->sample[i][j][k]);
                 hFront[j]->Fill(k, event->sample[i][j][k]);                 
               }
               break;
             case 1://Assuming 1 is back - to be checked
               for (int k = 0; k < NBSTRIPS; k++) { // strips   
                 hBackEnergy[j][k]->Fill(event->sample[i][j][k]);
                 hBack[j]->Fill(k, event->sample[i][j][k]);
               }
           }
         }
       }
     }
   }


   cout << "\r Processed all events" << endl;   

   // write on disk output file and close it
   outFile->Write();
   outFile->Close();
}
