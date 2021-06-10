/*****************************************************************************
 * Copyright (C) 2009-2021    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: A. Matta contact address: matta@lpccaen.in2p3.fr         *
 *                                                                           *
 * Creation Date  : May 2021                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  S034 analysis project                               *
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
#include"NPParticle.h"

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){

  Plastic= (TBigRIPSPlasticPhysics*) m_DetectorManager->GetDetector("BigRIPSPlastic");
  Minos= (TMinosPhysics*) m_DetectorManager->GetDetector("Minos");
  Nebula = (TNebulaPhysics*) m_DetectorManager->GetDetector("NEBULA");
  BDC = (TSamuraiBDCPhysics*) m_DetectorManager->GetDetector("SAMURAIBDC");
  FDC0 = (TSamuraiFDC0Physics*) m_DetectorManager->GetDetector("SAMURAIFDC0");
  FDC2 = (TSamuraiFDC2Physics*) m_DetectorManager->GetDetector("SAMURAIFDC2");
  Hodo = (TSamuraiHodoscopePhysics*) m_DetectorManager->GetDetector("SAMURAIHOD");

  FragmentTarget = NPL::EnergyLoss("He6_LH2.G4table","G4Table",1000 );
  m_field.LoadMap(30*deg,"field_map/180702-2,40T-3000.table.bin",10);
  m_field.SetFDC2Angle((59.930-90.0)*deg);
  m_field.SetFDC2R(FDC2->GetOffset().Z());
  InitOutputBranch();
  InitInputBranch();
  // for fdc/bdc alignement
  //file.open("Calibration/Pos/bdc.txt");
  He4=NPL::Particle("4He");
  He6=NPL::Particle("6He");
  He8=NPL::Particle("8He");
  He8.SetBeta(0.5152);
  n  =NPL::Particle("neutron");
  mhe6 = He6.Mass();
  mn   = n.Mass();
  sumM= mhe6+mn;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  Clear();
  
  static TLorentzVector LVHe6;
  static TLorentzVector LVn;
  Trigger=Trigger&0x00ff;


  // BigRIPS plastic
  double TF5=-1;
  double TF13=-1;
  double TF7=-1;
  double TOF_F5F13=-1;
  double TOF_F7F13=-1;
  unsigned int sizeP = Plastic->FP.size();
  for(unsigned int i = 0 ; i < sizeP ; i++){
    if(Plastic->FP[i]==5){
      TF5=Plastic->TSlew[i];
    }
    else if(Plastic->FP[i]==13&&Plastic->ID[i]==4){// two plastic at F13 taking only one
      TF13=Plastic->TSlew[i];
    }
    else if(Plastic->FP[i]==7){
      TF7=Plastic->TSlew[i];
    }
  }

  if(TF7>0 && TF13>0){
    // offset is adjusted to give the expected beta
    TOF_F7F13=TF13-TF7+6.71626e+02;
    static double LengthF7F13 = 117915-66409;
    Beta_b=(LengthF7F13/TOF_F7F13)/NPUNITS::c_light;
    // to find offset:
    //Beta_b=TOF_F7F13-LengthF7F13/He8.GetVelocity();
  }
  // Samurai-Minos
  if( Beta_b>0.5140  && Beta_b < 0.5165 // Correct Beta
      && Hodo->Charge.size()==1 && Hodo->Charge[0]>28 && Hodo->Charge[0]<42 && Hodo->Time[0]>58 && Hodo->Time[0]<68 // 6He in Hodo->cope
      && FDC2->PosX>-1500 && FDC2->PosX<1000 
      && FDC2->PosY>-500 && FDC2->PosY<500 
      && FDC0->PosX>-80 && FDC0->PosX<80 
      && FDC0->PosY>-80 && FDC0->PosY<80 // both FDC ok
      && Minos->Tracks_P0.size()==1     ) {  // p,pn only
    // Compute ThetaX and PhiY using Minos vertex and FDC0 X
    // Check if both BDC are reconstructed
    TVector3 BDC1=BDC->GetPos(1);
    TVector3 BDC2=BDC->GetPos(2);
    if( BDC1.Z()!=-10000 && BDC2.Z()!=-10000){
      TVector3 Vertex,delta;
      TVector3 P1 = Minos->Tracks_P0[0]+Minos->Tracks_Dir[0];

      // Vertex in Samurai Reference frame 
      MinimumDistanceTwoLines(BDC1,BDC2, 
          Minos->Tracks_P0[0], P1,
          Vertex, delta) ;
     
      TVector3 FDC0_Dir= FDC0->GetPos()-Vertex;
      FDC0_Dir=FDC0_Dir.Unit();
      TVector3 BDCDir=BDC2-BDC1;
      BDCDir=BDCDir.Unit();
      BDCDir*=(Vertex.Z()-BDC2.Z())/BDCDir.Z();
      BDCX=(BDC2+BDCDir).X();
      BDCY=(BDC2+BDCDir).Y();

      // XYZ relative to Minos entrance
      X=Vertex.X();
      Y=Vertex.Y();
      Z=Vertex.Z()+4650;

      // Fragment analysis 
      if(FDC0_Dir.Z()>0.6 && Z>0 && Z<150 && sqrt(X*X+Y*Y)<15){
        double FDC0_ThetaX = atan((FDC0->PosX-Vertex.X())/(1254.39-Z));
        double FDC0_PhiY   = atan((FDC0->PosY-Vertex.Y())/(1254.39-Z));
        Brho=m_field.FindBrho(FDC0->GetPos(),FDC0_Dir,FDC2->GetPos(),TVector3(0,0,1));
        He6.SetBrho(Brho);
        double Energy = He6.GetEnergy();
        if(Energy){
          Energy=FragmentTarget.EvaluateInitialEnergy(Energy,150-Z,FDC0_Dir.Angle(TVector3(0,0,1)));
          He6.SetKineticEnergy(Energy);
          Beta_f=He6.GetBeta();
          LVHe6.SetVectM(TVector3(0,0,0),mhe6);
          LVHe6.Boost(Beta_f*FDC0_Dir.Unit());
        }
        
      }
      
      // Neutron analysis
      if(Nebula->Charge.size()>0 && !Nebula->HasVeto()){
        // Index of first neutron hit
        unsigned int first = Nebula->GetFirstHit();
        TVector3 Pfirst = (Nebula->GetPos(first)-Vertex);
        double L = Pfirst.Mag();
        double TSBT= (Vertex.Z()+7377.56)/(Beta_b*c_light);
        TOF_n = Nebula->TOF[first]-TSBT-TF13;
        Beta_n = (L/TOF_n)/NPUNITS::c_light;
        LVn.SetVectM(TVector3(0,0,0),mn);
        LVn.Boost(Beta_n*Pfirst.Unit());
      }
      
      // a full event
      if(Beta_n&&Beta_f){
        He6.SetEnergyImpulsion(LVHe6);
        n.SetEnergyImpulsion(LVn);
        Erel = (LVHe6+LVn).Mag()-sumM;
      }
      // Calib//////////////////////////////////////////////////////////////////
 /*     static int count=0;
      if(Minos->Delta_Vertex < 5 && FDC2->PosX-252.55>0&&FDC0->GetPos().X()>-10000 && FDC0->Dir.Z()>0.9 && Minos->Z_Vertex>0&& sqrt(Minos->X_Vertex*Minos->X_Vertex+Minos->Y_Vertex*Minos->Y_Vertex)<15){
        file << FDC0->GetPos().X()   <<" " << FDC0->GetPos().Y() << " " << FDC0->GetPos().Z() <<" " ;
        file << Minos->X_Vertex      <<" " << Minos->Y_Vertex    << " " << Minos->Z_Vertex    << " " ;
        file << FDC2->GetPos().X()   <<" " << FDC2->GetPos().Y() << " " << FDC2->GetPos().Z() <<" " << FDC2->Dir.X() <<" " << FDC2->Dir.Y() << " " << FDC2->Dir.Z()<< endl;
        count ++;
      }
      if(count>1000)
        exit(1);
        */
    /*  static int count=0;
      if(Minos->Delta_Vertex < 5 && sqrt(Minos->X_Vertex*Minos->X_Vertex+Minos->Y_Vertex*Minos->Y_Vertex)<15 && Minos->Z_Vertex>-4650){
        file << BDC1.X()  <<" " << BDC1.Y()<< " " << BDC1.Z() <<" " ;
        file << BDC2.X()  <<" " << BDC2.Y()<< " " << BDC2.Z() <<" " ;
        file << Minos->X_Vertex      <<" " << Minos->Y_Vertex    << " " << Minos->Z_Vertex    << endl ;
        count ++;
      }
      if(count>10000)
        exit(1);
      */  

    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Clear(){
  Brho=-1000;
  BDCX=-1000;
  BDCY=-1000;
  Beta_f=-1000;
  Beta_n=-1000;
  Beta_b=-1000;
  X=Y=Z=-1000;
  TOF_n=-1000;
  Erel=-1000;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("Brho",&Brho,"Brho/D");
  RootOutput::getInstance()->GetTree()->Branch("BDCX",&BDCX,"BDCX/D");
  RootOutput::getInstance()->GetTree()->Branch("BDCY",&BDCY,"BDCY/D");
  RootOutput::getInstance()->GetTree()->Branch("X",&X,"X/D");
  RootOutput::getInstance()->GetTree()->Branch("Y",&Y,"Y/D");
  RootOutput::getInstance()->GetTree()->Branch("Z",&Z,"Z/D");
  RootOutput::getInstance()->GetTree()->Branch("Erel",&Erel,"Erel/D");
  RootOutput::getInstance()->GetTree()->Branch("Beta_f",&Beta_f,"Beta_f/D");
  RootOutput::getInstance()->GetTree()->Branch("Beta_n",&Beta_n,"Beta_n/D");
  RootOutput::getInstance()->GetTree()->Branch("Beta_b",&Beta_b,"Beta_b/D");
  RootOutput::getInstance()->GetTree()->Branch("TOF_n",&TOF_n,"TOF_n/D");
  RootOutput::getInstance()->GetTree()->Branch("Trigger",&Trigger,"Trigger/I");
} 
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  RootInput::getInstance()->GetChain()->SetBranchAddress("Trigger",&Trigger);
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
