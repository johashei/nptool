/***************************************************************************
 *   Charlie Paxman (cp00474@surrey.ac.uk)                                 *
 ***************************************************************************/

//-------------------------------
//C++
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
//Root
//#include <TVector3.h>
//NPTool
//#include "NPEnergyLoss.h"
//#include "NPReaction.h"
//#include "NPSystemOfUnits.h"
using namespace std;
//-------------------------------

void EventReader_AgataPhi(){

  // Read ROOT file and pull the tree	
  auto DataFile = new TFile("../../../../Outputs/Analysis/47Kdp_17Feb_AGATA-testing_WriteForMin_Large.root", "READ");
  auto PhysTree = (TTree*) DataFile->FindObjectAny("PhysicsTree");

  // Initilise the Mugast branch
  auto Mugast = new TMugastPhysics();

  // Initilise access variables for PhysTree
  static double T_MUGAST_VAMOS;
  static vector<double> 
	  *AGATA_GammaPx, 
	  *AGATA_GammaPy, 
	  *AGATA_GammaPz, 
	  *AGATA_GammaE, 
	  *AGATA_OrigBetaX, 
	  *AGATA_OrigBetaY, 
	  *AGATA_OrigBetaZ, 
	  *MugastHit;

  // Pull PhysTree branches
  auto MugVam_Branch = PhysTree->GetBranch("T_MUGAST_VAMOS");
  auto MugastHit_Branch = PhysTree->GetBranch("RawEnergy");
  auto AGATA_GammaPx_Branch = PhysTree->GetBranch("AGATA_GammaPx");
  auto AGATA_GammaPy_Branch = PhysTree->GetBranch("AGATA_GammaPy");
  auto AGATA_GammaPz_Branch = PhysTree->GetBranch("AGATA_GammaPz");
  auto AGATA_GammaE_Branch = PhysTree->GetBranch("AGATA_GammaE");
  auto AGATA_OrigBetaX_Branch = PhysTree->GetBranch("AGATA_OrigBetaX");
  auto AGATA_OrigBetaY_Branch = PhysTree->GetBranch("AGATA_OrigBetaY");
  auto AGATA_OrigBetaZ_Branch = PhysTree->GetBranch("AGATA_OrigBetaZ");

  // Set Mugast branch address
  PhysTree->SetBranchAddress("Mugast",&Mugast);

  // Set PhysTree variable addresses
  AGATA_GammaPx_Branch->SetAddress(&AGATA_GammaPx);
  AGATA_GammaPy_Branch->SetAddress(&AGATA_GammaPy);
  AGATA_GammaPz_Branch->SetAddress(&AGATA_GammaPz);
  AGATA_GammaE_Branch->SetAddress(&AGATA_GammaE);
  AGATA_OrigBetaX_Branch->SetAddress(&AGATA_OrigBetaX);
  AGATA_OrigBetaY_Branch->SetAddress(&AGATA_OrigBetaY);
  AGATA_OrigBetaZ_Branch->SetAddress(&AGATA_OrigBetaZ);
  MugVam_Branch->SetAddress(&T_MUGAST_VAMOS);
  MugastHit_Branch->SetAddress(&MugastHit);

  // Build loop variables
  unsigned int numEntries = PhysTree->GetEntries();
  unsigned int multiplicity = 0;

  // Open output file
  ofstream outfile;
  outfile.open("./RotBeta_17Feb_Large.txt");

  // Loop on entries
  for(unsigned int i=0; i<numEntries; i++){
    PhysTree->GetEntry(i);
    multiplicity = Mugast->TelescopeNumber.size(); 

    // Loop on MUGAST multiplicity
    //for(int m=0; m<multiplicity; m++){

//    cout << T_MUGAST_VAMOS << endl;
      // Gate on Timing
      if(abs(T_MUGAST_VAMOS-2777)<600){


        int gammaMultip = AGATA_GammaE->size();
        if(gammaMultip>=1){

          for(int g=0; g<gammaMultip; g++){

/*
            cout << AGATA_GammaPx->at(g)   << " " << endl;
            cout        << AGATA_GammaPy->at(g)   << " " << endl;
        cout        << AGATA_GammaPz->at(g)   << " " << endl;
        cout        << AGATA_GammaE->at(g)    << " " << endl;
        cout        << AGATA_OrigBetaX->at(g) << " " << endl;
        cout        << AGATA_OrigBetaY->at(g) << " " << endl;
        cout        << AGATA_OrigBetaZ->at(g) << endl; 
*/

            outfile << AGATA_GammaPx->at(g)   << " " 
                    << AGATA_GammaPy->at(g)   << " " 
                    << AGATA_GammaPz->at(g)   << " " 
                    << AGATA_GammaE->at(g)    << " " 
                    << AGATA_OrigBetaX->at(g) << " " 
                    << AGATA_OrigBetaY->at(g) << " " 
                    << AGATA_OrigBetaZ->at(g) << endl; 
	  }// for gamma mult
	}//if gamma mult

        //int gammaMultip = AddBack_EDC->size();
	  // gamma stuff (fails otherwise??)
          //if(gammaMultip!=0){
          //if(gammaMultip>=1){

	    // Gate on E_Gamma
	    //if(abs(AddBack_EDC->at(0)-0.143) < 0.002){
//	    if(abs(AddBack_EDC->at(0)-0.14) < 0.05){ // Only take gammas in range 0.09-0.19
//	      outfile << AddBack_EDC->at(0) << " " <<
//	                 AgataPhi->at(0) << " " <<
//	                 endl;
            //}//if 0.143
          //}//if gamma
      }//if timing
    //}//for m
  }//for i
  outfile.close();

  cout << " --- COMPLETE --- " << endl;
}
