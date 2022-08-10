/*****************************************************************************
 * Copyright (C) 2009-2014    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 * Edited for e793s                                                          *
 *                                                                           *
 * Creation Date  :                                                          *
 * Last update    : May 2021                                                 *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Class describing the property of an Analysis object                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include<iostream>
#include<string>

#include<fstream>

using namespace std;

#include "Analysis.h"
#include "NPFunction.h"
#include "NPAnalysisFactory.h"
#include "NPDetectorManager.h"
#include "NPOptionManager.h"
#include "macro/DefineColours.h"

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init() {
 ///////////////////////////////////////////////////////////////////////////////  

  if(NPOptionManager::getInstance()->HasDefinition("Sim")){
    cout << " == == == == SIMULATION == == == ==" << endl;
      isSim=true;
      isPhaseSpace=false;
  } else if (NPOptionManager::getInstance()->HasDefinition("Exp")) {
    cout << " == == == == EXPERIMENT == == == ==" << endl;
      isSim=false;
      isPhaseSpace=false;
  } else {
    cout << " == == == == PHASE SPACE == == == ==" << endl;
      isSim=false;
      //isSim=true;
      isPhaseSpace=true;
  }

  agata_zShift=51*mm;
  //BrhoRef=0.65;

  if(isSim && !isPhaseSpace){
//cout << "here_InIsSimLoop" << endl;
    Initial = new TInitialConditions();
    ReactionConditions = new TReactionConditions(); 
    RootInput::getInstance()->GetChain()->SetBranchAddress("InitialConditions",&Initial);

    // Initilize MGX histograms, because must be inside function
    string base1 = "ThetaCM_detected_MG";
    string base2 = "ThetaLab_detected_MG";
    for(int i=0; i<6; i++){
      int j=i+1;
      if (i==6){int j=7;}
      string name1 = base1 + to_string(j);	  
      string name2 = base2 + to_string(j);	  
      ThetaCM_detected_MGX[i] = new TH1F(name1.c_str(),name1.c_str(),NumThetaAngleBins,0,180); //900 bins for 0.2 angular bin width
      ThetaLab_detected_MGX[i] = new TH1F(name2.c_str(),name2.c_str(),NumThetaAngleBins,0,180);
    }

    // Initilize MMX histograms
    string base3 = "ThetaCM_detected_MM";
    string base4 = "ThetaLab_detected_MM";
    for(int i=0; i<5; i++){
      int j=i+1;
      string name1 = base3 + to_string(j);	  
      string name2 = base4 + to_string(j);	  
      ThetaCM_detected_MMX[i] = new TH1F(name1.c_str(),name1.c_str(),NumThetaAngleBins,0,180);
      ThetaLab_detected_MMX[i] = new TH1F(name2.c_str(),name2.c_str(),NumThetaAngleBins,0,180);
    }
  }

  // initialize input and output branches
  InitOutputBranch();
  InitInputBranch();
  
  M2 = (TMust2Physics*)  m_DetectorManager -> GetDetector("M2Telescope");
  MG = (TMugastPhysics*) m_DetectorManager -> GetDetector("Mugast");
  ML = (TModularLeafPhysics*) m_DetectorManager -> GetDetector("ModularLeaf"); 
//  CATS = (TCATSPhysics*) m_DetectorManager->GetDetector("CATSDetector");

  // get reaction information
  reaction.ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  OriginalBeamEnergy = reaction.GetBeamEnergy();
  reaction.Print(); //TESTING PARSER

  // get beam position from .reaction file
  Beam = (NPL::Beam*) reaction.GetParticle1(); 
  XBeam = Beam->GetMeanX();
  YBeam = Beam->GetMeanY();

  cout << " ---------------- Beam Position ---------------- " << endl;
  cout << " \tX = " << XBeam << " mm\tY = " << YBeam  << " mm" << endl;
  cout << " ----------------------------------------------- " << endl;

  // target thickness
  TargetThickness = m_DetectorManager->GetTargetThickness();
  string TargetMaterial = m_DetectorManager->GetTargetMaterial();
  
  cout << " ---------------- Target Values ---------------- " << endl;
  cout << " \tZ = " << m_DetectorManager->GetTargetZ() << " mm\tT = " << TargetThickness  << " mm" << endl;
  cout << " ----------------------------------------------- " << endl;
  
  // energy losses
  string light = NPL::ChangeNameToG4Standard(reaction.GetParticle3()->GetName());
  string beam = NPL::ChangeNameToG4Standard(reaction.GetParticle1()->GetName());
  LightTarget = NPL::EnergyLoss(light+"_"+TargetMaterial+".G4table","G4Table",100 );
  LightAl = NPL::EnergyLoss(light+"_Al.G4table" ,"G4Table",100);
  LightSi = NPL::EnergyLoss(light+"_Si.G4table" ,"G4Table",100);
  BeamTargetELoss = NPL::EnergyLoss(beam+"_"+TargetMaterial+".G4table","G4Table",100);

  FinalBeamEnergy = BeamTargetELoss.Slow(OriginalBeamEnergy, 0.5*TargetThickness, 0);
  reaction.SetBeamEnergy(FinalBeamEnergy); 

  cout << "\033[91m Running for reaction "
       << reaction.GetParticle1()->GetName() << "("
       << reaction.GetParticle2()->GetName() << ","
       << reaction.GetParticle3()->GetName() << ")"
       << reaction.GetParticle4()->GetName() << endl;

  cout << "\033[36m Beam energy at mid-target: " << FinalBeamEnergy << "\033[37m"<< endl;

  // initialize various parameters
  Rand = TRandom3();
  ParticleMult = 0;
  GammaMult = 0;
  DetectorNumber = 0;
  ThetaNormalTarget = 0;
  ThetaM2Surface = 0;
  ThetaMGSurface = 0;
  Si_E_M2 = 0;
  CsI_E_M2 = 0;
  Energy = 0;
  ThetaGDSurface = 0;
  elab_tmp = 0;
  thetalab_tmp = 0;
  philab_tmp = 0;
  ThetaGDSurface = 0;
  X.clear();
  Y.clear();
  Z.clear();
  dE = 0;
  BeamDirection = TVector3(0,0,1);
  //nbTrack=0;
  //nbHits=0;
  //count=0;
  AHeavy=reaction.GetParticle4()->GetA();
  ALight=reaction.GetParticle3()->GetA(); 
  MHeavy=reaction.GetParticle4()->Mass();
  MLight=reaction.GetParticle3()->Mass();
  //bool writetoscreen=true;

  for(int i=0;i<GATCONF_SIZE;i++){ // loop over the bits
    GATCONF_Counter[i] = 0 ; 
  }

//  ThetaCM_detected->Sumw2();
//  ThetaLab_detected->Sumw2();
//cout << "here_endInit" << endl;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  // Reinitiate calculated variable
//cout << "here_TreatEvent" << endl;
  ReInitValue();
  if(isSim && !isPhaseSpace){
    ThetaCM_emmitted->Fill(ReactionConditions->GetThetaCM());
    ThetaLab_emmitted->Fill(ReactionConditions->GetTheta(0));
  }

  GATCONF_MASTER=ML->GetCalibratedValue("GATCONF_MASTER");
  for(int i=0;i<GATCONF_SIZE;i++){ // loop over the bits
    if(GATCONF_MASTER & (unsigned int)pow(2,i)){ // test if ith bit is on
      GATCONF_Counter[i]++ ; // increment the array
      //      cout <<  std::bitset<16>(GATCONF_MASTER) << " " << i+1 << endl;
    }
  }

  if(isSim && !isPhaseSpace){
    OriginalELab = ReactionConditions->GetKineticEnergy(0);
    OriginalThetaLab = ReactionConditions->GetTheta(0);
    BeamEnergy = ReactionConditions->GetBeamEnergy();
  }

  TVector3 BeamDirection(0.,0.,1.);
  BeamImpact = TVector3(XBeam,YBeam,m_DetectorManager->GetTargetZ()); 

  //ParticleMult=M2->Si_E.size();////+MG->DSSD_E.size();
  ParticleMult=M2->Si_E.size()+MG->DSSD_E.size();
//cout << "here_BeforeMustLoop" << endl;
  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////////// LOOP on MUST2  ////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  unsigned int sizeM2 = M2->Si_E.size();
  for(unsigned int countMust2 = 0; countMust2 < sizeM2; countMust2++){
//cout << "here_InMustLoop" << endl;

    /************************************************/
    // Part 0 : Get the useful Data
    // MUST2
    int TelescopeNumber = M2->TelescopeNumber[countMust2];

    if(isSim && !isPhaseSpace){
      ThetaCM_detected_MM->Fill(ReactionConditions->GetThetaCM());
      ThetaLab_detected_MM->Fill(ReactionConditions->GetTheta(0));

      int MMX = TelescopeNumber-1;
      ThetaCM_detected_MMX[MMX]->Fill(ReactionConditions->GetThetaCM());
      ThetaLab_detected_MMX[MMX]->Fill(ReactionConditions->GetTheta(0));
    }

    /************************************************/
    // Part 1 : Impact Angle
    ThetaM2Surface = 0;
    ThetaNormalTarget = 0;
    thetalab_tmp = 0;
    philab_tmp = 0;
    TVector3 HitDirection = M2 -> GetPositionOfInteraction(countMust2) - BeamImpact ;
    thetalab_tmp = HitDirection.Angle(BeamDirection);
    philab_tmp = HitDirection.Phi();

    X.push_back(M2->GetPositionOfInteraction(countMust2).X());
    Y.push_back(M2->GetPositionOfInteraction(countMust2).Y());
    Z.push_back(M2->GetPositionOfInteraction(countMust2).Z());

    ThetaM2Surface = HitDirection.Angle(- M2->GetTelescopeNormal(countMust2) );
    ThetaNormalTarget = HitDirection.Angle( TVector3(0.0,0.0,1) ) ;

    /************************************************/
    // Part 2 : Impact Energy
    Energy = 0;
    elab_tmp = 0;
    Si_E_M2 = M2->Si_E[countMust2];
    CsI_E_M2= M2->CsI_E[countMust2];

    // if CsI
    if(CsI_E_M2>0 ){
      // The energy in CsI is calculate form dE/dx Table because
      Energy = CsI_E_M2;
      //if(!isSim){
      Energy = LightAl.EvaluateInitialEnergy(Energy, 0.4*micrometer, ThetaM2Surface);
      //}
      Energy+=Si_E_M2;
    }
    else
      Energy = Si_E_M2;

    RawEnergy.push_back(Energy);

//    if(!isSim){
      // Evaluate energy using the thickness
      elab_tmp = LightAl.EvaluateInitialEnergy(
		      Energy, 
		      0.4*micrometer, 
		      ThetaM2Surface);
      ELoss_Al.push_back(Energy-elab_tmp);
      double elab_tmp2 = elab_tmp;
      // Target Correction
      elab_tmp = LightTarget.EvaluateInitialEnergy(
		      elab_tmp, 
		      0.5*TargetThickness, 
		      ThetaNormalTarget);
      ELoss_Target.push_back(elab_tmp2-elab_tmp);
      ELoss.push_back(Energy-elab_tmp);
//    } else {elab_tmp = Energy;}

    ELab.push_back(elab_tmp);

    /************************************************/

    /************************************************/
    // Part 3 : Excitation Energy Calculation
    Ex.push_back(reaction.ReconstructRelativistic(elab_tmp,thetalab_tmp));
    Ecm.push_back(Energy*(AHeavy+ALight)/(4*AHeavy*cos(thetalab_tmp)*cos(thetalab_tmp)));
    /************************************************/

    /************************************************/
    // Part 4 : Theta CM Calculation
    
    ThetaCM.push_back(reaction.EnergyLabToThetaCM(elab_tmp, thetalab_tmp)/deg);
    /************************************************/

    ThetaLab.push_back(thetalab_tmp/deg);
    PhiLab.push_back(philab_tmp/deg);
//cout << "here_EndMustLoop" << endl;
  }

//cout << "here_BeforeMugastLoop" << endl;
  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////////// LOOP on MUGAST ////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  unsigned int sizeMG = MG->DSSD_E.size();
  for(unsigned int countMugast = 0; countMugast<sizeMG; countMugast++){


    if(isSim && !isPhaseSpace){
      int MGX = MG->TelescopeNumber[0];
      MGX = MGX-1;
      if(MGX==6){MGX=5;}

      ThetaCM_detected_MG->Fill(ReactionConditions->GetThetaCM());
      ThetaCM_detected_MGX[MGX]->Fill(ReactionConditions->GetThetaCM());

      ThetaLab_detected_MG->Fill(ReactionConditions->GetTheta(0));
      ThetaLab_detected_MGX[MGX]->Fill(ReactionConditions->GetTheta(0));
    }

    // Part 1 : Impact Angle
    ThetaMGSurface = 0;
    ThetaNormalTarget = 0;
    thetalab_tmp = 0;
    philab_tmp = 0;
    TVector3 HitDirection = MG->GetPositionOfInteraction(countMugast) - BeamImpact;
    //TVector3 tempVector;
    //if(MG->TelescopeNumber[0]==3){
    //  tempVector = {MG->GetPositionOfInteraction(countMugast).X(),
    //                MG->GetPositionOfInteraction(countMugast).Y(),
    //                MG->GetPositionOfInteraction(countMugast).Z()+1.0}; //add 1mm to MG3
    //  HitDirection = tempVector - BeamImpact;
    //  if(warning){cout << "!!! EDITING MG3 !!!" << endl; warning=false;}
    //}
    thetalab_tmp = HitDirection.Angle(BeamDirection);
    philab_tmp = HitDirection.Phi();

    X.push_back( MG -> GetPositionOfInteraction(countMugast).X());
    Y.push_back( MG -> GetPositionOfInteraction(countMugast).Y());
    Z.push_back( MG -> GetPositionOfInteraction(countMugast).Z());

    ThetaMGSurface = HitDirection.Angle( MG -> GetTelescopeNormal(countMugast) );
    ThetaNormalTarget = HitDirection.Angle( TVector3(0.0,0.0,1.0) ) ;

    // Part 2 : Impact Energy
    Energy = elab_tmp = 0;
    Energy = MG->GetEnergyDeposit(countMugast);
    RawEnergy.push_back(Energy);

    if(!isSim){
      elab_tmp = LightAl.EvaluateInitialEnergy(
		    Energy,              //particle energy after Al
		    0.4*micrometer,      //thickness of Al
		    ThetaMGSurface);     //angle of impingement
      ELoss_Al.push_back(Energy-elab_tmp);
      double elab_tmp2 = elab_tmp;
      elab_tmp = LightTarget.EvaluateInitialEnergy(
		    elab_tmp,            //particle energy after leaving target
		    TargetThickness*0.5, //distance passed through target
		    ThetaNormalTarget);  //angle of exit from target
      ELoss_Target.push_back(elab_tmp2-elab_tmp);
      ELoss.push_back(Energy-elab_tmp);
    } else { //TESTING DIFFERENT ENERGY LOSSES IN SIMULATION
      elab_tmp = Energy; //so I can add and remove sections
      //elab_tmp = LightSi.EvaluateInitialEnergy(
      //	    elab_tmp,            //particle energy after Si
      //	    0.5*500.*micrometer, //thickness of Si
      //	    ThetaMGSurface);     //angle of impingement
      //elab_tmp = LightAl.EvaluateInitialEnergy(
      //              elab_tmp,            //particle energy after Al
      //	            0.4*micrometer,      //thickness of Al
      //              ThetaMGSurface);     //angle of impingement
      elab_tmp = LightTarget.EvaluateInitialEnergy(
    		    elab_tmp,            //particle energy after leaving target
    		    TargetThickness*0.5, //distance passed through target
    		    ThetaNormalTarget);  //angle of exit from target
    }

    ELab.push_back(elab_tmp);

    //cout << "===============" << endl;
    //cout << "RawE:\t" << RawEnergy.back() << endl;
    //cout << "ELAl:\t" << ELoss_Al.back() << endl;
    //cout << "ELCD:\t" << ELoss_Target.back() << endl;
    //cout << "ELTt:\t" << ELoss.back() << endl;
    //cout << "ELab:\t" << ELab.back() << endl;


    // Part 3 : Excitation Energy Calculation
    //if(!isSim){ //TESTING!!!!
      Ex.push_back(reaction.ReconstructRelativistic(elab_tmp,thetalab_tmp));
      Ecm.push_back(elab_tmp*(AHeavy+ALight)/(4*AHeavy*cos(thetalab_tmp)*cos(thetalab_tmp)));
    //}

    // Part 4 : Theta CM Calculation
    ThetaLab.push_back(thetalab_tmp/deg);
    PhiLab.push_back(philab_tmp/deg);
    ThetaCM.push_back(reaction.EnergyLabToThetaCM(elab_tmp, thetalab_tmp)/deg);

    if(sizeMG==1){
      MG_T = MG->DSSD_T[0];
      MG_E = MG->DSSD_E[0];
      MG_X = MG->DSSD_X[0];
      MG_Y = MG->DSSD_Y[0];
      MG_D = MG->TelescopeNumber[0];
    }

  }//end loop Mugast




  ////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////// LOOP on AGATA ////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  // Agata by track
  /*
  for(int j=0; j<nbTrack; j++){ // all multiplicity
    TLorentzVector GammaLV;
    // Measured E
    double Egamma=trackE[j]/1000.; // From keV to MeV 
    // Gamma detection position
    // TrackZ1 to be corrected there is a shift of +51mm
    TVector3 GammaHit(trackX1[j],trackY1[j],trackZ1[j]+agata_zShift); 
    // TVector3 GammaHit(trackX1[0],trackY1[0],trackZ1[0]); 
    // Gamma Direction 
    TVector3 GammaDirection = GammaHit-BeamImpact;
    GammaDirection = GammaDirection.Unit();
    // Beta from Two body kinematic
    //TVector3 beta = reaction.GetEnergyImpulsionLab_4().BoostVector();
    // Beta from the Beam mid target 
    reaction.GetKinematicLine4();
    TVector3 beta(0,0,-reaction.GetNucleus4()->GetBeta());

    double ThetaGamma = GammaDirection.Angle(BeamDirection)/deg;
    // Construct LV in lab frame
    GammaLV.SetPx(Egamma*GammaDirection.X());
    GammaLV.SetPy(Egamma*GammaDirection.Y());
    GammaLV.SetPz(Egamma*GammaDirection.Z());
    GammaLV.SetE(Egamma);
    // Boost back in CM
    GammaLV.Boost(beta);
    // Get EDC
    EDC.push_back(GammaLV.Energy());
  }
  */

  // Agata add back is not always multiplicity 1 ?? NO, not necessarily!
  for(int i=0; i<nbAdd; i++){
  //if(nbAdd==1){
    TLorentzVector GammaLV;
    // Measured E
    double Egamma = AddE[i]/1000.; // From keV to MeV 
    // Gamma detection position
    // TrackZ1 to be corrected there is a shift of +51mm
    TVector3 GammaHit(AddX[i],AddY[i],AddZ[i]+agata_zShift); 
    // TVector3 GammaHit(trackX1[0],trackY1[0],trackZ1[0]); 


//    cout << "gDir " << GammaHit.X() << " "
//                    << GammaHit.Y() << " "
//                    << GammaHit.Z() << " ";
    /* TEST adding constant to Phi */
    //double phitemp = GammaHit.Phi();
    //double phiconst = -2.0;
    //GammaHit.SetPhi(phitemp+phiconst);
//    GammaHit.RotateZ(0.0);

//    cout << "gDir " << GammaHit.X() << " "
//                    << GammaHit.Y() << " "
//                    << GammaHit.Z() << endl;

    // Gamma Direction 
    TVector3 GammaDirection = GammaHit-BeamImpact;
    GammaDirection = GammaDirection.Unit();


    /* Beta from Two body kinematic */
    TVector3 beta = reaction.GetEnergyImpulsionLab_4().BoostVector();
    //TVector3 beta = reaction->GetEnergyImpulsionLab_4().BoostVector();

//    cout << "bDir " << beta.X() << " "
//                    << beta.Y() << " "
//                    << beta.Z() ;
//    beta.RotateX(0.002847); beta.RotateY(3.144869); beta.RotateZ(0.095923);
    beta.RotateY(M_PI);
//    cout << "bDir " << beta.X() << " "
//                    << beta.Y() << " "
//                    << beta.Z() << endl;

    /* Original beta */
//    TVector3 beta(0,0,-0.1257);

    // For beta rotation minimization
    AGATA_GammaPx.push_back(Egamma*GammaDirection.X());
    AGATA_GammaPy.push_back(Egamma*GammaDirection.Y());
    AGATA_GammaPz.push_back(Egamma*GammaDirection.Z());
    AGATA_GammaE.push_back(Egamma);
    AGATA_OrigBetaX.push_back(beta.X());
    AGATA_OrigBetaY.push_back(beta.Y());
    AGATA_OrigBetaZ.push_back(beta.Z());

    /* Other fills */
    double ThetaGamma = GammaDirection.Angle(BeamDirection)/deg;
    // Construct LV in lab frame
    GammaLV.SetPx(Egamma*GammaDirection.X());
    GammaLV.SetPy(Egamma*GammaDirection.Y());
    GammaLV.SetPz(Egamma*GammaDirection.Z());
    GammaLV.SetE(Egamma);
    // Boost back in CM
    GammaLV.Boost(beta);
    // Get EDC
    AddBack_EDC.push_back(GammaLV.Energy());
    AddBack_EDC2.push_back(GammaLV.Energy());

  }

  float offset[21] = {0.,27.4,46.7,41.1,71.2,58.7,32.1,40.4,-46.3,15.1,9.9,-39.1,-39.8,-5.3,18.5,0.7,-28.5,-7,-22,-8.5,-8.5};
  float slope[21] = {1.,0.96,0.93,0.93,0.89,0.91,0.95,0.93,1.06,0.97,0.99,1.05,1.04,1.01,0.97,1.,1.03,0.99,1.,1.02,1.02};

  for(unsigned int i=0; i<MW_Nr; i++){
    MWT[i] = -1500;
    for(unsigned short j=0; j<21; j++)
      if(MW_N[i]==j)
        MWT[i] = (MW_T[i]-offset[j])/slope[j];
  }


}

////////////////////////////////////////////////////////////////////////////////
void FillSolidAngles(TH1F* hSA, TH1F* hDet, TH1F* hEmm){
  if(!filledCline){
    cout << RED << "FILLING CLINE VECTOR! SHOULD ONLY OCCUR ONCE!" << RESET << endl;
    for(int t=0; t<NumThetaAngleBins; t++){
      double angleMin = (t)*(180./NumThetaAngleBins);
      double angleMax = (t+1)*(180./NumThetaAngleBins);
      cout << " Angle " << angleMin 
	      << " to " << angleMax
	      << endl; 
      clineVal.push_back(2.0*M_PI*(cos(angleMin*degtorad) - cos(angleMax*degtorad)));
      clineX.push_back(angleMin+((angleMax-angleMin)/2.0));
      filledCline=true;
    }
  }

  for (int i = 0; i < hSA->GetNbinsX(); ++i){
    double val = hDet->GetBinContent(i) / hEmm->GetBinContent(i);
    double valerr = val * sqrt( 
      pow(hDet->GetBinError(i) / hDet->GetBinContent(i), 2) +
      pow(hEmm->GetBinError(i) / hEmm->GetBinContent(i), 2) );
    if (isnan(val)) { val = 0; valerr = 0; }
    val *= clineVal.at(i);
    valerr *= clineVal.at(i);

    hSA->SetBinContent(i, val);
    hSA->SetBinError(i, valerr);
  }
}


//void DivideByCline(TH1F* histo){
//
//  for(int b=0; b<NumThetaAngleBins; b++){
//    int store = histo->GetBinContent(b);
//    histo->SetBinContent(b,(int) store/clineVal.at(b));
///  }
//}



////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
  cout << endl << "\e[1;32m GATCONF Statistics " << endl ;
  for(int i=0;i<GATCONF_SIZE;i++){ // loop over the bits
    cout << GATCONF_Label[i]  << "\t\t" << GATCONF_Counter[i] << endl ; 
  }
  cout << endl ;

  if(isSim && !isPhaseSpace){
    
    //TObjArray HistList(0);
    TList *HistList = new TList();

    // MUST 2 DETECTOR #5
    auto Efficiency_CM_MM = new TH1F(*ThetaCM_detected_MM);
    Efficiency_CM_MM->SetName("Efficiency_CM_MM");
    Efficiency_CM_MM->SetTitle("Efficiency_CM_MM");
    Efficiency_CM_MM->Sumw2();
    auto Efficiency_Lab_MM = new TH1F(*ThetaLab_detected_MM);
    Efficiency_Lab_MM->SetName("Efficiency_Lab_MM");
    Efficiency_Lab_MM->SetTitle("Efficiency_Lab_MM");
    Efficiency_Lab_MM->Sumw2();
     
    auto SolidAngle_CM_MM = new TH1F(*ThetaCM_detected_MM);
    SolidAngle_CM_MM->SetName("SolidAngle_CM_MM");
    SolidAngle_CM_MM->SetTitle("SolidAngle_CM_MM");
    auto SolidAngle_Lab_MM = new TH1F(*ThetaLab_detected_MM);
    SolidAngle_Lab_MM->SetName("SolidAngle_Lab_MM");
    SolidAngle_Lab_MM->SetTitle("SolidAngle_Lab_MM");

    Efficiency_CM_MM->Divide(ThetaCM_emmitted);
    Efficiency_Lab_MM->Divide(ThetaLab_emmitted);

    //double dt_MM = 180./Efficiency_Lab_MM->GetNbinsX();
    //cout << "Angular infinitesimal (MM) = " << dt_MM << "deg " << endl;
    //auto Cline_MM = new TF1("Cline_MM",Form("1./(2*%f*sin(x*%f/180.)*%f*%f/180.)",M_PI,M_PI,dt_MM,M_PI),0,180);

    FillSolidAngles(SolidAngle_CM_MM, ThetaCM_detected_MM, ThetaCM_emmitted);
    FillSolidAngles(SolidAngle_Lab_MM, ThetaLab_detected_MM, ThetaLab_emmitted);
    
    //SolidAngle_CM_MM->Divide(ThetaCM_emmitted);
    //SolidAngle_CM_MM->Divide(Cline_MM,1);
    //SolidAngle_CM_MM->Divide(Cline);
    //SolidAngle_Lab_MM->Divide(ThetaLab_emmitted);
    //SolidAngle_Lab_MM->Divide(Cline_MM,1);
    //SolidAngle_Lab_MM->Divide(Cline);

    HistList->Add(ThetaCM_emmitted);
    HistList->Add(ThetaLab_emmitted);
    HistList->Add(ThetaCM_detected_MM);
    HistList->Add(ThetaLab_detected_MM);
    HistList->Add(Efficiency_CM_MM);
    HistList->Add(Efficiency_Lab_MM);
    HistList->Add(SolidAngle_CM_MM);
    HistList->Add(SolidAngle_Lab_MM);
    //HistList->Add(Cline_MM);
    //HistList->Add(Cline);


    // MUGAST
    auto Efficiency_CM_MG = new TH1F(*ThetaCM_detected_MG);
    Efficiency_CM_MG->SetName("Efficiency_CM_MG");
    Efficiency_CM_MG->SetTitle("Efficiency_CM_MG");
    Efficiency_CM_MG->Sumw2();
    auto Efficiency_Lab_MG = new TH1F(*ThetaLab_detected_MG);
    Efficiency_Lab_MG->SetName("Efficiency_Lab_MG");
    Efficiency_Lab_MG->SetTitle("Efficiency_Lab_MG");
    Efficiency_Lab_MG->Sumw2();
     
    auto SolidAngle_CM_MG = new TH1F(*ThetaCM_detected_MG);
    SolidAngle_CM_MG->SetName("SolidAngle_CM_MG");
    SolidAngle_CM_MG->SetTitle("SolidAngle_CM_MG");
    auto SolidAngle_Lab_MG = new TH1F(*ThetaLab_detected_MG);
    SolidAngle_Lab_MG->SetName("SolidAngle_Lab_MG");
    SolidAngle_Lab_MG->SetTitle("SolidAngle_Lab_MG");

    Efficiency_CM_MG->Divide(ThetaCM_emmitted);
    Efficiency_Lab_MG->Divide(ThetaLab_emmitted);

    //double dt_MG = 180./Efficiency_Lab_MG->GetNbinsX();
    //cout << "Angular infinitesimal (MG) = " << dt_MG << "deg " << endl;
    //auto Cline_MG = new TF1("Cline_MG",Form("1./(2*%f*sin(x*%f/180.)*%f*%f/180.)",M_PI,M_PI,dt_MG,M_PI),0,180);

    /* Testing method for better errors in SolidAngle histograms */
    FillSolidAngles(SolidAngle_CM_MG, ThetaCM_detected_MG, ThetaCM_emmitted);
    //SolidAngle_CM_MG->Divide(Cline_MG,1);
    //SolidAngle_CM_MG->Divide(Cline);
    FillSolidAngles(SolidAngle_Lab_MG, ThetaLab_detected_MG, ThetaLab_emmitted);
    //SolidAngle_Lab_MG->Divide(Cline_MG,1);
    //SolidAngle_Lab_MG->Divide(Cline);

    HistList->Add(ThetaCM_detected_MG);
    HistList->Add(ThetaLab_detected_MG);
    HistList->Add(Efficiency_CM_MG);
    HistList->Add(Efficiency_Lab_MG);
    HistList->Add(SolidAngle_CM_MG);
    HistList->Add(SolidAngle_Lab_MG);
    //HistList->Add(Cline_MG);

    // MUGAST INDIVIDUAL
    TH1F *SolidAngle_CM_MGX[6];
    string base1 = "SolidAngle_CM_MG";
    for(int i=0;i<6;i++){
      int j=i+1;
      if(j==6){int j=7;}
      string name = base1 + to_string(j);	  
      SolidAngle_CM_MGX[i] = new TH1F(*ThetaCM_detected_MGX[i]);
      SolidAngle_CM_MGX[i]->SetName(name.c_str());
      SolidAngle_CM_MGX[i]->SetTitle(name.c_str());
      FillSolidAngles(SolidAngle_CM_MGX[i], ThetaCM_detected_MGX[i], ThetaCM_emmitted);
      //SolidAngle_CM_MGX[i]->Divide(Cline_MG,1);
      //SolidAngle_CM_MGX[i]->Divide(Cline);
    }

    TH1F *SolidAngle_Lab_MGX[6];
    string base2 = "SolidAngle_Lab_MG";
    for(int i=0;i<6;i++){
      int j=i+1;
      if(j==6){int j=7;}
      string name = base2 + to_string(j);	  
      SolidAngle_Lab_MGX[i] = new TH1F(*ThetaLab_detected_MGX[i]);
      SolidAngle_Lab_MGX[i]->SetName(name.c_str());
      SolidAngle_Lab_MGX[i]->SetTitle(name.c_str());
      FillSolidAngles(SolidAngle_Lab_MGX[i], ThetaLab_detected_MGX[i], ThetaLab_emmitted);
      //SolidAngle_Lab_MGX[i]->Divide(Cline_MG,1);
      //SolidAngle_Lab_MGX[i]->Divide(Cline);
    }

    // MUST2 INDIVIDUAL
    TH1F *SolidAngle_CM_MMX[6];
    string base3 = "SolidAngle_CM_MM";
    for(int i=0;i<5;i++){
      int j=i+1;
      string name = base3 + to_string(j);	  
      SolidAngle_CM_MMX[i] = new TH1F(*ThetaCM_detected_MMX[i]);
      SolidAngle_CM_MMX[i]->SetName(name.c_str());
      SolidAngle_CM_MMX[i]->SetTitle(name.c_str());
      FillSolidAngles(SolidAngle_CM_MMX[i], ThetaCM_detected_MMX[i], ThetaCM_emmitted);
      //SolidAngle_CM_MMX[i]->Divide(Cline_MM,1);
      //SolidAngle_CM_MMX[i]->Divide(Cline);
    }

    TH1F *SolidAngle_Lab_MMX[6];
    string base4 = "SolidAngle_Lab_MM";
    for(int i=0;i<5;i++){
      int j=i+1;
      string name = base4 + to_string(j);	  
      SolidAngle_Lab_MMX[i] = new TH1F(*ThetaLab_detected_MMX[i]);
      SolidAngle_Lab_MMX[i]->SetName(name.c_str());
      SolidAngle_Lab_MMX[i]->SetTitle(name.c_str());
      FillSolidAngles(SolidAngle_Lab_MMX[i], ThetaLab_detected_MMX[i], ThetaLab_emmitted);
      //SolidAngle_Lab_MMX[i]->Divide(Cline_MM,1);
      //SolidAngle_Lab_MMX[i]->Divide(Cline);
    }

    for(int i=0; i<6; i++){HistList->Add(ThetaCM_detected_MGX[i]);}
    for(int i=0; i<6; i++){HistList->Add(ThetaLab_detected_MGX[i]);}
    for(int i=0; i<5; i++){HistList->Add(ThetaCM_detected_MMX[i]);}
    for(int i=0; i<5; i++){HistList->Add(ThetaLab_detected_MMX[i]);}

    for(int i=0; i<6; i++){HistList->Add(SolidAngle_CM_MGX[i]);}
    for(int i=0; i<6; i++){HistList->Add(SolidAngle_Lab_MGX[i]);}
    for(int i=0; i<5; i++){HistList->Add(SolidAngle_CM_MMX[i]);}
    for(int i=0; i<5; i++){HistList->Add(SolidAngle_Lab_MMX[i]);}
    
    auto clineValGraph = new TGraph(NumThetaAngleBins);
    clineValGraph->SetName("clineValGraph");
    clineValGraph->SetTitle("clineValGraph");
    for(int b=0; b<NumThetaAngleBins; b++){
      clineValGraph->SetPoint(b,clineX.at(b),clineVal.at(b));
    }
    HistList->Add(clineValGraph);

    auto HistoFile = new TFile("~/Programs/nptool/Projects/e793s/SolidAngle_HistFile_New.root","RECREATE");
    HistList->Write();
    HistoFile->Close();


    cout << "!!! MAKE SURE YOU RENAME THE SOLID ANGLE FILE! !!!" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
//cout << "here_InitOutput" << endl;
  //RootOutput::getInstance()->GetTree()->Branch("Ex",&Ex,"Ex/D");
  RootOutput::getInstance()->GetTree()->Branch("Ex",&Ex);
  //RootOutput::getInstance()->GetTree()->Branch("EDC",&EDC,"EDC/D");
  RootOutput::getInstance()->GetTree()->Branch("EDC",&EDC);
  RootOutput::getInstance()->GetTree()->Branch("AddBack_EDC",&AddBack_EDC);
  RootOutput::getInstance()->GetTree()->Branch("AddBack_EDC2",&AddBack_EDC2);
  RootOutput::getInstance()->GetTree()->Branch("AGATA_GammaPx",&AGATA_GammaPx);
  RootOutput::getInstance()->GetTree()->Branch("AGATA_GammaPy",&AGATA_GammaPy);
  RootOutput::getInstance()->GetTree()->Branch("AGATA_GammaPz",&AGATA_GammaPz);
  RootOutput::getInstance()->GetTree()->Branch("AGATA_GammaE",&AGATA_GammaE);
  RootOutput::getInstance()->GetTree()->Branch("AGATA_OrigBetaX",&AGATA_OrigBetaX);
  RootOutput::getInstance()->GetTree()->Branch("AGATA_OrigBetaY",&AGATA_OrigBetaY);
  RootOutput::getInstance()->GetTree()->Branch("AGATA_OrigBetaZ",&AGATA_OrigBetaZ);
  RootOutput::getInstance()->GetTree()->Branch("EAgata",&EAgata,"EAgata/D");
  RootOutput::getInstance()->GetTree()->Branch("ELab",&ELab);
  RootOutput::getInstance()->GetTree()->Branch("Ecm",&Ecm);
  RootOutput::getInstance()->GetTree()->Branch("RawEnergy",&RawEnergy);
  RootOutput::getInstance()->GetTree()->Branch("ELoss_Al",&ELoss_Al);
  RootOutput::getInstance()->GetTree()->Branch("ELoss_Target",&ELoss_Target);
  RootOutput::getInstance()->GetTree()->Branch("ELoss",&ELoss);
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab);
  RootOutput::getInstance()->GetTree()->Branch("PhiLab",&PhiLab);
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM",&ThetaCM);
  RootOutput::getInstance()->GetTree()->Branch("Run",&Run,"Run/I");
  RootOutput::getInstance()->GetTree()->Branch("X",&X);
  RootOutput::getInstance()->GetTree()->Branch("Y",&Y);
  RootOutput::getInstance()->GetTree()->Branch("Z",&Z);
  RootOutput::getInstance()->GetTree()->Branch("dE",&dE,"dE/D");
  RootOutput::getInstance()->GetTree()->Branch("MG_T",MG_T);
  RootOutput::getInstance()->GetTree()->Branch("MG_E",MG_E);
  RootOutput::getInstance()->GetTree()->Branch("MG_X",MG_X);
  RootOutput::getInstance()->GetTree()->Branch("MG_Y",MG_Y);
  RootOutput::getInstance()->GetTree()->Branch("MG_D",MG_D);

  // Vamos 
  RootOutput::getInstance()->GetTree()->Branch("LTS",&LTS,"LTS/l");
/*    RootOutput::getInstance()->GetTree()->Branch("T_FPMW_CATS1",&T_FPMW_CATS1,"T_FPMW_CATS1/s");
    RootOutput::getInstance()->GetTree()->Branch("T_FPMW_CATS2_C",&T_FPMW_CATS2_C,"T_FPMW_CATS2_C/F");
    RootOutput::getInstance()->GetTree()->Branch("T_FPMW_HF",&T_FPMW_HF,"T_FPMW_CATS1/s");
    RootOutput::getInstance()->GetTree()->Branch("T_FPMW_HF_C",&T_FPMW_HF_C,"T_FPMW_CATS1/F");
    RootOutput::getInstance()->GetTree()->Branch ("IC",IC, "IC[12]/F");
    RootOutput::getInstance()->GetTree()->Branch ("ICRawPUMult01",ICRawPUMult01, "ICRawPUMult01[7]/s");

    RootOutput::getInstance()->GetTree()->Branch("T_CATS2_HF",&T_CATS2_HF,"T_CATS2_HF/s");
    RootOutput::getInstance()->GetTree()->Branch("T_MUGAST_FPMW",&T_MUGAST_FPMW,"T_MUGAST_FPMW/s");
    RootOutput::getInstance()->GetTree()->Branch("T_FPMW_CATS2",&T_FPMW_CATS2,"T_FPMW_CATS2/s");
    RootOutput::getInstance()->GetTree()->Branch("T_FPMW_HF",&T_FPMW_HF,"T_FPMW_HF/s");

    // DC0
    RootOutput::getInstance()->GetTree()->Branch ("DC0_X", &DC0_X, "DC0_X/F");
    RootOutput::getInstance()->GetTree()->Branch ("DC0_Y", &DC0_Y, "DC0_Y/F");
    // DC1
    RootOutput::getInstance()->GetTree()->Branch ("DC1_X", &DC1_X, "DC1_X/F");
    RootOutput::getInstance()->GetTree()->Branch ("DC1_Y", &DC1_Y, "DC1_Y/F");
    // DC2
    RootOutput::getInstance()->GetTree()->Branch ("DC2_X", &DC2_X, "DC2_X/F");
    RootOutput::getInstance()->GetTree()->Branch ("DC2_Y", &DC2_Y, "DC2_Y/F");
    // DC3
    RootOutput::getInstance()->GetTree()->Branch ("DC3_X", &DC3_X, "DC3_X/F");
    RootOutput::getInstance()->GetTree()->Branch ("DC3_Y", &DC3_Y, "DC3_Y/F");

    RootOutput::getInstance()->GetTree()->Branch ("Xf", &Xf, "Xf/F");
    RootOutput::getInstance()->GetTree()->Branch ("Tf", &Tf, "Tf/F");

    RootOutput::getInstance()->GetTree()->Branch ("Brho", &Brho, "Brho/F");
    RootOutput::getInstance()->GetTree()->Branch ("Theta", &Theta, "Theta/F");
    RootOutput::getInstance()->GetTree()->Branch ("Phi", &Phi, "Phi/F");
    RootOutput::getInstance()->GetTree()->Branch ("Path", &Path, "Path/F");

    RootOutput::getInstance()->GetTree()->Branch ("EWIRE_1_1", &EWIRE_1_1, "EWIRE_1_1/s");
    RootOutput::getInstance()->GetTree()->Branch ("EWIRE_1_2", &EWIRE_1_2, "EWIRE_1_2/s");
    RootOutput::getInstance()->GetTree()->Branch ("EWIRE_2_1", &EWIRE_2_1, "EWIRE_2_1/s");
    RootOutput::getInstance()->GetTree()->Branch ("EWIRE_2_2", &EWIRE_2_2, "EWIRE_2_2/s");
*/
  RootOutput::getInstance()->GetTree()->Branch ("MW_Nr", &MW_Nr, "MW_Nr/I");
  RootOutput::getInstance()->GetTree()->Branch ("MW_T", MW_T, "MW_T[MW_Nr]/F");
  RootOutput::getInstance()->GetTree()->Branch ("MW_N", MW_N, "MW_N[MW_Nr]/s");
  RootOutput::getInstance()->GetTree()->Branch ("MWT", MWT, "MWT[MW_Nr]/F");

  RootOutput::getInstance()->GetTree()->Branch ("AGAVA_VAMOSTS", &AGAVA_VAMOSTS, "AGAVA_VAMOSTS/l");
/*    RootOutput::getInstance()->GetTree()->Branch("mT",&mT,"mT/D");
    RootOutput::getInstance()->GetTree()->Branch("mV",&mV,"mV/D");
    RootOutput::getInstance()->GetTree()->Branch("mD",&mD,"mD/D");
    RootOutput::getInstance()->GetTree()->Branch("mBeta",&mBeta,"mBeta/D");
    RootOutput::getInstance()->GetTree()->Branch("mGamma",&mGamma,"mGamma/D");
    RootOutput::getInstance()->GetTree()->Branch("mM_Q",&mM_Q,"mM_Q/D");
    RootOutput::getInstance()->GetTree()->Branch("mM",&mM,"mM/D");
    RootOutput::getInstance()->GetTree()->Branch("mE",&mE,"mE/D");
    RootOutput::getInstance()->GetTree()->Branch("mdE",&mdE,"mdE/D");
*/

    // DC0->3_X/_Y Xf, Tf, Brho, Path, EWIRE_1_1 all, MW_Nr, MW_T, AGAVA_VAMOSTS, 
    // Agata
    // Time stamp of the agata trigger
    /*
    RootOutput::getInstance()->GetTree()->Branch("TStrack",&TStrack,"TStrack/l");

    // Array of reconstructed tracks
    RootOutput::getInstance()->GetTree()->Branch("nbTrack",&nbTrack,"nbTrack/I");
    RootOutput::getInstance()->GetTree()->Branch("trackE",trackE,"trackE[nbTrack]/F");
    RootOutput::getInstance()->GetTree()->Branch("trackX1",trackX1,"trackX1[nbTrack]/F");
    RootOutput::getInstance()->GetTree()->Branch("trackY1",trackY1,"trackY1[nbTrack]/F");
    RootOutput::getInstance()->GetTree()->Branch("trackZ1",trackZ1,"trackZ1[nbTrack]/F");
    RootOutput::getInstance()->GetTree()->Branch("trackT",trackT,"trackT[nbTrack]/F");
    RootOutput::getInstance()->GetTree()->Branch("trackCrystalID",trackCrystalID,"trackCrystalID[nbTrack]/I");
    */
    
  //Agata Addback
  RootOutput::getInstance()->GetTree()->Branch("nbAdd",&nbAdd,"nbAdd/I");
  RootOutput::getInstance()->GetTree()->Branch("TSHit",&TSHit,"TSHit/l");
  RootOutput::getInstance()->GetTree()->Branch("AddE",AddE,"AddE[nbAdd]/F");
  RootOutput::getInstance()->GetTree()->Branch("AddX",AddX,"AddX[nbAdd]/F");
  RootOutput::getInstance()->GetTree()->Branch("AddY",AddY,"AddY[nbAdd]/F");
  RootOutput::getInstance()->GetTree()->Branch("AddZ",AddZ,"AddZ[nbAdd]/F");
  
  RootOutput::getInstance()->GetTree()->Branch("EventWeight",&EventWeight,"EventWeight/D");

  if(isSim && !isPhaseSpace){
    RootOutput::getInstance()->GetTree()->Branch("OriginalELab",
	      &OriginalELab,"OriginalELab/D");
    RootOutput::getInstance()->GetTree()->Branch("OriginalThetaLab",
	      &OriginalThetaLab,"OriginalThetaLab/D");
    RootOutput::getInstance()->GetTree()->Branch("BeamEnergy",
	      &BeamEnergy,"BeamEnergy/D");
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){

//cout << "here_InitInput" << endl;
  SetBranchStatus();
  // RootInput:: getInstance()->GetChain()->SetBranchAddress("GATCONF",&vGATCONF);
  //
  // Vamos
  RootInput::getInstance()->GetChain()->SetBranchAddress("LTS",&LTS);
/*    RootInput::getInstance()->GetChain()->SetBranchAddress("T_FPMW_CATS1",&T_FPMW_CATS1);
    RootInput::getInstance()->GetChain()->SetBranchAddress("T_FPMW_CATS2_C",&T_FPMW_CATS2_C);
    RootInput::getInstance()->GetChain()->SetBranchAddress("T_FPMW_HF_C",&T_FPMW_HF_C);

    RootInput::getInstance()->GetChain()->SetBranchAddress("T_CATS2_HF",&T_CATS2_HF);
    RootInput::getInstance()->GetChain()->SetBranchAddress("T_MUGAST_FPMW",&T_MUGAST_FPMW);
    RootInput::getInstance()->GetChain()->SetBranchAddress("T_FPMW_CATS2",&T_FPMW_CATS2);
    RootInput::getInstance()->GetChain()->SetBranchAddress("T_FPMW_HF",&T_FPMW_HF);

    RootInput::getInstance()->GetChain()->SetBranchAddress("Xf", &Xf);
    RootInput::getInstance()->GetChain()->SetBranchAddress("Tf", &Tf);

    RootInput::getInstance()->GetChain()->SetBranchAddress("Brho", &Brho);
    RootInput::getInstance()->GetChain()->SetBranchAddress("Theta", &Theta);
    RootInput::getInstance()->GetChain()->SetBranchAddress("Phi", &Phi);
    RootInput::getInstance()->GetChain()->SetBranchAddress("Path", &Path);

    RootInput::getInstance()->GetChain()->SetBranchAddress("EWIRE_1_1", &EWIRE_1_1);
    RootInput::getInstance()->GetChain()->SetBranchAddress("EWIRE_1_2", &EWIRE_1_2);
    RootInput::getInstance()->GetChain()->SetBranchAddress("EWIRE_2_1", &EWIRE_2_1);
    RootInput::getInstance()->GetChain()->SetBranchAddress("EWIRE_2_2", &EWIRE_2_2);
*/
  RootInput::getInstance()->GetChain()->SetBranchAddress("MWTVM", &MW_Nr);
  RootInput::getInstance()->GetChain()->SetBranchAddress("MWTV", MW_T);
  RootInput::getInstance()->GetChain()->SetBranchAddress("MWTVN", MW_N);
  RootInput::getInstance()->GetChain()->SetBranchAddress("AGAVA_VAMOSTS", &AGAVA_VAMOSTS);

  // Agata
     /*
    RootInput::getInstance()->GetChain()->SetBranchAddress("TStrack",&TStrack);
    RootInput::getInstance()->GetChain()->SetBranchAddress("nbTrack",&nbTrack);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackE",trackE);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackX1",trackX1);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackY1",trackY1);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackZ1",trackZ1);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackT",trackT);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackCrystalID",trackCrystalID);
    */
  RootInput::getInstance()->GetChain()->SetBranchAddress("TSHit",&TSHit);
  RootInput::getInstance()->GetChain()->SetBranchAddress("nbAdd",&nbAdd);
  RootInput::getInstance()->GetChain()->SetBranchAddress("AddE",AddE);
  RootInput::getInstance()->GetChain()->SetBranchAddress("AddX",AddX);
  RootInput::getInstance()->GetChain()->SetBranchAddress("AddY",AddY);
  RootInput::getInstance()->GetChain()->SetBranchAddress("AddZ",AddZ);

  if(isPhaseSpace){
    RootInput:: getInstance()->GetChain()->SetBranchAddress("EventWeight",&EventWeight);
  }

  if(isSim && !isPhaseSpace){
    //RootInput:: getInstance()->GetChain()->SetBranchStatus("InitialConditions",true );
    //RootInput:: getInstance()->GetChain()->SetBranchStatus("fIC_*",true );
    RootInput::getInstance()->GetChain()->SetBranchAddress("InitialConditions",
		    &Initial);
    RootInput::getInstance()->GetChain()->SetBranchAddress("InteractionCoordinates",
		    &Interaction);
    //RootInput:: getInstance()->GetChain()->SetBranchStatus("ReactionConditions",true );
    //RootInput:: getInstance()->GetChain()->SetBranchStatus("fRC_*",true );
    RootInput::getInstance()->GetChain()->SetBranchAddress("ReactionConditions",
		    &ReactionConditions);
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::SetBranchStatus(){
//cout << "here_SetBranchStatus" << endl;
  // Set Branch status 
  RootInput::getInstance()->GetChain()->SetBranchStatus("LTS",true);
 /*   RootInput::getInstance()->GetChain()->SetBranchStatus("T_FPMW_CATS1",true);
    RootInput::getInstance()->GetChain()->SetBranchStatus("T_FPMW_CATS2_C",true);
    RootInput::getInstance()->GetChain()->SetBranchStatus("T_FPMW_HF",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("T_FPMW_HF_C",true);
    
    RootInput::getInstance()->GetChain()->SetBranchStatus("T_CATS2_HF",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("T_MUGAST_FPMW",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("T_FPMW_CATS2",true );

    RootInput::getInstance()->GetChain()->SetBranchStatus("Xf",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("Tf",true );

    RootInput::getInstance()->GetChain()->SetBranchStatus("Brho",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("Path",true );

    RootInput::getInstance()->GetChain()->SetBranchStatus("EWIRE_1_1",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("EWIRE_1_2",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("EWIRE_2_1",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("EWIRE_2_2",true );
*/
  RootInput::getInstance()->GetChain()->SetBranchStatus("MW_Nr",true );
  RootInput::getInstance()->GetChain()->SetBranchStatus("MW_T",true  );

  RootInput::getInstance()->GetChain()->SetBranchStatus("MWTV",true );
  RootInput::getInstance()->GetChain()->SetBranchStatus("MWTVN",true  );
  RootInput::getInstance()->GetChain()->SetBranchStatus("MWTVM",true  );

  RootInput::getInstance()->GetChain()->SetBranchStatus("AGAVA_VAMOSTS",true  );
  // Agata
    /*
    RootInput::getInstance()->GetChain()->SetBranchStatus("TStrack",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("nbTrack",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("trackE",true  );
    RootInput::getInstance()->GetChain()->SetBranchStatus("trackX1",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("trackY1",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("trackZ1",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("trackT",true  );
    RootInput::getInstance()->GetChain()->SetBranchStatus("trackCrystalID",true);
    */
  RootInput::getInstance()->GetChain()->SetBranchStatus("TSHit",true );
  RootInput::getInstance()->GetChain()->SetBranchStatus("nbAdd",true );
  RootInput::getInstance()->GetChain()->SetBranchStatus("AddE",true  );
  RootInput::getInstance()->GetChain()->SetBranchStatus("AddX",true  );
  RootInput::getInstance()->GetChain()->SetBranchStatus("AddY",true  );
  RootInput::getInstance()->GetChain()->SetBranchStatus("AddZ",true  );

  if(isPhaseSpace){
    RootInput:: getInstance()->GetChain()->SetBranchStatus("EventWeight",true );
  }

  if(isSim && !isPhaseSpace){
    RootInput:: getInstance()->GetChain()->SetBranchStatus("InitialConditions",true );
    RootInput:: getInstance()->GetChain()->SetBranchStatus("fIC_*",true );
    RootInput:: getInstance()->GetChain()->SetBranchStatus("InteractionCoordinates",true );
    RootInput:: getInstance()->GetChain()->SetBranchStatus("ReactionConditions",true );
    RootInput:: getInstance()->GetChain()->SetBranchStatus("fRC_*",true );
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
//cout << "here_ReInit" << endl;
  Ex.clear();
  Ecm.clear();
  AddBack_EDC.clear();
  AddBack_EDC2.clear();
  AddBack_EDC2.push_back(-1.0); //offset by 1
  AGATA_GammaPx.clear();
  AGATA_GammaPy.clear();
  AGATA_GammaPz.clear();
  AGATA_GammaE.clear();
  AGATA_OrigBetaX.clear();
  AGATA_OrigBetaY.clear();
  AGATA_OrigBetaZ.clear();
  EAgata = -1000;
  ELab.clear();
  RawEnergy.clear(); 
  ELoss_Al.clear();
  ELoss_Target.clear();
  ELoss.clear();
  ThetaLab.clear();
  PhiLab.clear();
  ThetaCM.clear();
  //relative_angle = -1000;
  //sum_angle = -1000;
  X.clear();
  Y.clear();
  Z.clear();
  dE= -1000;
  ParticleMult = 0;
  GammaMult = 0;
  DetectorNumber = 0;
  ThetaNormalTarget = 0;
  ThetaM2Surface = 0;
  ThetaMGSurface = 0;
  Si_E_M2 = 0;
  CsI_E_M2 = 0;
  Energy = 0;
  ThetaGDSurface = 0;
  BeamDirection = TVector3(0,0,1);
  elab_tmp = 0;
  thetalab_tmp = 0;
  philab_tmp = 0;

  filledCline=false;

  MG_T=-1000;
  MG_E=-1000;
  MG_X=-1000;
  MG_Y=-1000;
  MG_D=-1000;

  //EDC= -1000;
  EDC.clear();
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the AnalysisFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct(){
  return (NPL::VAnalysis*) new Analysis();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_analysis{
    public:
      proxy_analysis(){
        NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct);
      }
  };

  proxy_analysis p_analysis;
}
