/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: XAUTHORX  contact address: XMAILX                        *
 *                                                                           *
 * Creation Date  : XMONTHX XYEARX                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Scone analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include<iostream>
using namespace std;
#include"Analysis.h"
#include"NPAnalysisFactory.h"
#include"NPDetectorManager.h"
#include"RootInput.h"
#include"RootOutput.h"

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  Scone= (TSconePhysics*) m_DetectorManager->GetDetector("Scone");
  InitialConditions = new TInitialConditions();

  m_DetectedNeutron = 0;

  InitInputBranch();
  InitOutputBranch();

  m_entries = RootInput::getInstance()->GetChain()->GetEntries();

  E_init = 0;
  E_new = 0;
  new_energy = false;

  vE_init.clear();
  vDetectedNeutron.clear();
}


////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  ReInitValue();
  if(InitialConditions->GetParticleMultiplicity()>0){
    E_init = InitialConditions->GetKineticEnergy(0);
  }

  if(E_new == E_init) new_energy = false;
  if(E_new != E_init){
    E_new = E_init;
    vE_init.push_back(E_new);
    new_energy = true;
  }

  if(new_energy == true){
    if(m_DetectedNeutron!=0) vDetectedNeutron.push_back(m_DetectedNeutron);
    m_DetectedNeutron = 0;
  }

  E_sum = 0;
  if(Scone->Energy.size()>0){
    for(int i=0; i<Scone->Energy.size(); i++){
      if(Scone->Energy[i]>E_max) E_max = Scone->Energy[i];
      if(Scone->Time[i]>Time_max) Time_max = Scone->Time[i];
      E_sum += Scone->Energy[i];
    }
  }
  //E_sum = E_sum - E_init;
  //if(Time_max>50) m_DetectedNeutron++;
  if(Time_max>50 && E_sum>0) m_DetectedNeutron++;
  
  if(Scone->GammaEnergy.size()>0) E_sum_gamma = 0;
  else E_sum_gamma = -100;
  for(int i=0; i<Scone->GammaEnergy.size(); i++){
    E_sum_gamma += Scone->GammaEnergy[i];
  }
  
  if(Scone->ProtonEnergy.size()>0) E_sum_proton = 0;
  else E_sum_proton = -100;
  for(int i=0; i<Scone->ProtonEnergy.size(); i++){
    E_sum_proton += Scone->ProtonEnergy[i];
  }
  if(Scone->ProtonEnergy.size()>0) E_mean_proton = E_sum_proton/Scone->ProtonEnergy.size();

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  E_init = -10000;
  E_max = -10000;
  Time_max = -10000;
  E_mean_proton = -100;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch()
{
  RootOutput::getInstance()->GetTree()->Branch("E_init",&E_init,"E_init/D");
  RootOutput::getInstance()->GetTree()->Branch("E_max",&E_max,"E_max/D");
  RootOutput::getInstance()->GetTree()->Branch("E_sum",&E_sum,"E_sum/D");
  RootOutput::getInstance()->GetTree()->Branch("E_sum_gamma",&E_sum_gamma,"E_sum_gamma/D");
  RootOutput::getInstance()->GetTree()->Branch("E_sum_proton",&E_sum_proton,"E_sum_proton/D");
  RootOutput::getInstance()->GetTree()->Branch("E_mean_proton",&E_mean_proton,"E_mean_proton/D");
  RootOutput::getInstance()->GetTree()->Branch("Time_max",&Time_max,"Time_max/D");
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  RootInput::getInstance()->GetChain()->SetBranchStatus("InitialConditions",true);
  RootInput::getInstance()->GetChain()->SetBranchStatus("fIC_*",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("InitialConditions",&InitialConditions);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
  // For last energy //
  vDetectedNeutron.push_back(m_DetectedNeutron);
  
  /*cout << "Number of Init energy treated: " << vE_init.size() << endl;
  cout << "With initial energy: " << endl;
  for(int i=0; i< vE_init.size(); i++)
    cout << "* " << vE_init[i] << endl;

  cout << "Size of vDetectedNeutron: " << vDetectedNeutron.size() << endl;
  cout << "DetectedNeutron: " << endl;

  ofstream ofile;
  ofile.open("macro/eff_scone_test.txt");
  //ofile.open("macro/eff_scone_natGd25um.txt");
  for(int i=0; i< vDetectedNeutron.size(); i++){
    //cout << "* " << vE_init[i] << " / " << vDetectedNeutron[i]/vDetectedNeutron[0]*99.4 << endl;
    cout << "* " << vE_init[i] << " / " << vDetectedNeutron[i]/1e5*100 << endl;
    //ofile << vE_init[i] << "  " << vDetectedNeutron[i]/vDetectedNeutron[0]*99.4 << endl;
    ofile << vE_init[i] << "  " << vDetectedNeutron[i]/1e5*100 << endl;
  }
  ofile.close();*/
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct(){
  return (NPL::VAnalysis*) new Analysis();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy{
    public:
      proxy(){
        NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct);
      }
  };

  proxy p;
}

