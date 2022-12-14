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

void EventReader(){

  // Read ROOT file and pull the tree	
  //auto DataFile = new TFile("../../../../Outputs/Analysis/47K_Full_08July.root", "READ");
  auto DataFile = new TFile("../../../../Outputs/Analysis/47Kdp_11Apr22_PartII.root", "READ");
  auto PhysTree = (TTree*) DataFile->FindObjectAny("PhysicsTree");

  // Initilise the Mugast branch
  auto Mugast = new TMugastPhysics();

  // Initilise access variables for PhysTree
  static double T_MUGAST_VAMOS;
  static vector<double> *X, *Y, *Z, *RawEnergy, *AddBack_EDC, *ThetaLab;

  // Pull PhysTree branches
  auto Energy_Branch = PhysTree->GetBranch("RawEnergy");
  auto Gamma_Branch = PhysTree->GetBranch("AddBack_EDC");
  auto X_Branch = PhysTree->GetBranch("X");
  auto Y_Branch = PhysTree->GetBranch("Y");
  auto Z_Branch = PhysTree->GetBranch("Z");
  auto MugVam_Branch = PhysTree->GetBranch("T_MUGAST_VAMOS");
  auto ThetaLab_Branch = PhysTree->GetBranch("ThetaLab");

  // Set Mugast branch address
  PhysTree->SetBranchAddress("Mugast",&Mugast);

  // Set PhysTree variable addresses
  Energy_Branch->SetAddress(&RawEnergy);
  Gamma_Branch->SetAddress(&AddBack_EDC);
  X_Branch->SetAddress(&X);
  Y_Branch->SetAddress(&Y);
  Z_Branch->SetAddress(&Z);
  MugVam_Branch->SetAddress(&T_MUGAST_VAMOS);
  ThetaLab_Branch->SetAddress(&ThetaLab);

  // Build loop variables
  unsigned int numEntries = PhysTree->GetEntries();
  unsigned int multiplicity = 0;

  // Open output file
  ofstream outfile;
  //outfile.open("./XYZE_ShiftMG3_Sep08_GateOn0143.txt");
  outfile.open("./XYZE_GateOn1838_11Apr22_PtII.txt");

  // Loop on entries
  for(unsigned int i=0; i<numEntries; i++){
    PhysTree->GetEntry(i);
    multiplicity = Mugast->TelescopeNumber.size(); 

    // Loop on MUGAST multiplicity
    for(int m=0; m<multiplicity; m++){

      // Gate on Timing
      if(abs(T_MUGAST_VAMOS-2700)<400){
        int gammaMultip = AddBack_EDC->size();
//        cout << "GAMMA MULT " << gammaMultip << endl;

        //if(ThetaLab->at(m) <= 130. && ThetaLab->at(m) >= 100.){

	  // gamma stuff (fails otherwise??)
          //if(gammaMultip!=0){
          if(gammaMultip>=1){
            //for(int g; g<gammaMultip; g++){
	
	    // Gate on E_Gamma
	    if(abs(AddBack_EDC->at(0)-1.838) < 0.005){
	      outfile << Mugast->TelescopeNumber.at(m) << " " <<
	                 X->at(m) << " " <<
                         Y->at(m) << " " <<
	                 Z->at(m) << " " <<
	                 RawEnergy->at(m) << endl;
            }//if 0.143
          //}//for g
          }//if gamma
	//}//if theta (new!! testing!!!
      }//if timing
    }//for m
  }//for i
  outfile.close();

  cout << " --- COMPLETE --- " << endl;
}
