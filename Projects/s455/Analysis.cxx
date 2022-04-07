/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: P. Morfouace contact address: pierre.morfouace2@cea.fr   *
 *                                                                           *
 * Creation Date  : June 2021                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Sofia analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include<iostream>
#include<algorithm>
using namespace std;
#include"Analysis.h"
#include"NPAnalysisFactory.h"
#include"NPDetectorManager.h"
#include"NPPhysicalConstants.h"
#include"NPGlobalSystemOfUnits.h"


////////////////////////////////////////////////////////////////////////////////
struct TofPair
{
  double x = -1000;
  double y = -1000;
  double tof = -1000;
  double velocity = -1;
  double beta = -1;
  double theta_in = -10;
  double theta_out = -10;
  int plastic = -1;
  //
  int isLorR  = -1;
  // *** isLorR = 1 -> Left
  // *** isLorR = 2 -> Right 
  //
  int isUorD  = -1;
  // *** isUorD = 1 -> Up
  // *** isUorD = 2 -> Down
  int section = -1;
  double Esec = -1;
  double DT = -100;
  double x3 = -1000;
};


////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  SofBeamID = new TSofBeamID();
  SofFF = new TSofFissionFragment();
  SofSci= (TSofSciPhysics*) m_DetectorManager->GetDetector("SofSci");
  SofTrim= (TSofTrimPhysics*) m_DetectorManager->GetDetector("SofTrim");
  SofTwim= (TSofTwimPhysics*) m_DetectorManager->GetDetector("SofTwim");
  SofTofW= (TSofTofWPhysics*) m_DetectorManager->GetDetector("SofTofW");
  SofAt= (TSofAtPhysics*) m_DetectorManager->GetDetector("SofAt");
  SofMwpc= (TSofMwpcPhysics*) m_DetectorManager->GetDetector("SofMwpc");

  InitParameter();
  InitOutputBranch();
  LoadSpline();

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  ReInitValue();
  //cout << "************" << endl;
  RunID = fRunID;
  BeamAnalysis();

  unsigned int sofsci_size = SofSci->DetectorNbr.size();
  if(sofsci_size==2){
    double start_time = SofSci->TimeNs[1];
    SofTofW->SetTofAlignedValue(36);
    SofTofW->SetStartTime(start_time);
    SofTofW->BuildPhysicalEvent();

    FissionFragmentAnalysis();
    //BeamFragmentAnalysis();
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::BeamFragmentAnalysis(){
  unsigned int softofw_size = SofTofW->PlasticNbr.size();

  double L_CC = 8.45;
  TofPair TofHit;
  if(softofw_size==1){
    TofHit.plastic  = SofTofW->PlasticNbr[0];
    TofHit.x        = SofTofW->CalPosX[0];
    TofHit.y        = SofTofW->CalPosY[0];
    TofHit.tof      = SofTofW->CalTof[0];
    TofHit.velocity = L_CC/TofHit.tof;
    TofHit.beta     = TofHit.velocity * m/ns / NPUNITS::c_light;

    double Brho = 9.62543 + 0.0076642*TofHit.x;
    double Lfactor = 9.17/L_CC;
    double Beta = TofHit.beta*Lfactor;
    double Gamma1 = 1. / sqrt(1 - Beta * Beta);

    double AoQ = Brho / (3.10761 * Beta * Gamma1);

    SofFF->SetTOF(TofHit.tof);
    SofFF->SetTofPosX(TofHit.x);
    SofFF->SetTofPosY(TofHit.y);
    SofFF->SetPlastic(TofHit.plastic);

    SofFF->SetBeta(Beta);
    SofFF->SetGamma(Gamma1);
    SofFF->SetAoQ(AoQ);
    SofFF->SetBrho(Brho);


  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::FissionFragmentAnalysis(){
  unsigned int softofw_size = SofTofW->PlasticNbr.size();
  unsigned int softwim_size = SofTwim->SectionNbr.size();

  double E1 = -1;
  double E2 = -1;
  double E3 = -1;
  double E4 = -1;
  double DT1 = -1000;
  double DT2 = -1000;
  double DT3 = -1000;
  double DT4 = -1000;
  double Theta1 = -1000;
  double Theta2 = -1000;
  double Theta3 = -1000;
  double Theta4 = -1000;
  double L_CC = 8.45;
  double Beta_norm = 0.745;
  double Gamma1 = -1;
  double Gamma2 = -1;
  double Beta_Z1 = -1;
  double Beta_Z2 = -1;
  double Brho1 = -1;
  double Brho2 = -1;
  double AoQ1 = -1;
  double AoQ2 = -1;
  double A1 = -1;
  double A2 = -1;
  double Z1 = -1;
  double Z2 = -1;
  double Zsum = -1;
  int iZ1 = -1;
  int iZ2 = -1;
  int iZsum = -1;

  TofPair TofHit[2];
  if(softofw_size==2){
    for(unsigned int i=0; i<softofw_size; i++){
      TofHit[i].plastic  = SofTofW->PlasticNbr[i];
      TofHit[i].x        = SofTofW->CalPosX[i];
      TofHit[i].y        = SofTofW->CalPosY[i];
      TofHit[i].tof      = SofTofW->CalTof[i];
      TofHit[i].velocity = L_CC/TofHit[i].tof;
      TofHit[i].beta     = TofHit[i].velocity * m/ns / NPUNITS::c_light;
    }

    if(TofHit[0].x>TofHit[1].x){
      TofHit[0].isLorR = 1;
      TofHit[1].isLorR = 2;
    }
    else if(TofHit[0].x<TofHit[1].x){
      TofHit[0].isLorR = 2;
      TofHit[1].isLorR = 1;
    }

    if(TofHit[0].y>TofHit[1].y){
      TofHit[0].isUorD = 1;
      TofHit[1].isUorD = 2;
    }
    else if(TofHit[0].y<TofHit[1].y){
      TofHit[0].isUorD = 2;
      TofHit[1].isUorD = 1;
    }
  }

  vector<double> X1;
  vector<double> X2;
  vector<double> X3;
  vector<double> Y1;
  vector<double> Y2;
  vector<double> Y3;
  for(unsigned int i=0; i<SofMwpc->DetectorNbr.size(); i++){
    if(SofMwpc->DetectorNbr[i]==2){
      SofFF->SetPosX1(SofMwpc->PositionX1[i]);
      SofFF->SetPosY1(SofMwpc->PositionY[i]);
    }
    if(SofMwpc->DetectorNbr[i]==3){
      SofFF->SetPosX2(SofMwpc->PositionX1[i]);
      SofFF->SetPosY2(SofMwpc->PositionY[i]);
    }
    if(SofMwpc->DetectorNbr[i]==4){
      if(SofMwpc->PositionX1[i]!=-1000)
        X3.push_back(SofMwpc->PositionX1[i]);
      
      if(SofMwpc->PositionY[i]!=-1000)
        Y3.push_back(SofMwpc->PositionY[i]);

      //SofFF->SetPosX3(SofMwpc->PositionX1[i]);
      //SofFF->SetPosY3(SofMwpc->PositionY[i]);
    }
  }
  

  for(unsigned int i=0; i<2; i++){
    double tofx = TofHit[i].x;
    for(unsigned int k=0; k<X3.size(); k++){
      double posx = X3[k];
      if(abs(tofx-posx) < 150){
        if(abs(tofx-posx)<abs(tofx-TofHit[i].x3))
          TofHit[i].x3 = posx;
      }
    }
  }

  int mult1 = SofTwim->mult1;
  int mult2 = SofTwim->mult2;
  int mult3 = SofTwim->mult3;
  int mult4 = SofTwim->mult4;

  int multL = mult1 + mult2;
  int multR = mult3 + mult4;

  if(softwim_size>1){
    if( (mult1>1 && mult1<17) || (mult2>1 && mult2<17) || (mult3>1 && mult3<17) || (mult4>1 && mult4<17)){
      for(unsigned int i=0; i< softwim_size; i++){
        int sec = SofTwim->SectionNbr[i];

        if(sec==1){
          E1 = SofTwim->EnergySection[i];
          DT1 = SofTwim->DriftTime[i];
          Theta1 = SofTwim->Theta[i];
        }
        if(sec==2){
          E2 = SofTwim->EnergySection[i];
          DT2 = SofTwim->DriftTime[i];
          Theta2 = SofTwim->Theta[i];
        }
        if(sec==3){
          E3 = SofTwim->EnergySection[i]; 
          DT3 = SofTwim->DriftTime[i];
          Theta3 = SofTwim->Theta[i];
        }
        if(sec==4){
          E4 = SofTwim->EnergySection[i];     
          DT4 = SofTwim->DriftTime[i];
          Theta4 = SofTwim->Theta[i];
        }
      }

      if(softwim_size>2){
        if(E1>0 && E2>0){
          E1 = E1+E2;
          E2 = -1;
        }
        if(E3>0 && E4>0){
          E3 = E3+E4;
          E4 = -1;
        }
      }

      if(E1>0)
        E1 = E1/16;
      if(E2>0)
        E2 = E2/16;
      if(E3>0)
        E3 = E3/16;
      if(E4>0)
        E4 = E4/16;


      // *** case 1 *** //
      if(E1!=-1 && E2!=-1){
        if(TofHit[0].isUorD==1 && TofHit[1].isUorD==2){
          TofHit[0].Esec = E2;
          TofHit[1].Esec = E1;

          TofHit[0].theta_in = Theta2;
          TofHit[1].theta_in = Theta1;
 
          TofHit[0].DT = DT2;
          TofHit[1].DT = DT1;

          TofHit[0].section = 2;
          TofHit[1].section = 1;
        }
        if(TofHit[0].isUorD==2 && TofHit[1].isUorD==1){
          TofHit[0].Esec = E1;
          TofHit[1].Esec = E2;

          TofHit[0].theta_in = Theta1;
          TofHit[1].theta_in = Theta2;
 
          TofHit[0].DT = DT1;
          TofHit[1].DT = DT2;

          TofHit[0].section = 1;
          TofHit[1].section = 2;
        }
      }

      // *** case 2 *** //
      if(E1!=-1 && E3!=-1){
        if(TofHit[0].isUorD==1 && TofHit[1].isUorD==2){
          TofHit[0].Esec = E3;
          TofHit[1].Esec = E1;

          TofHit[0].theta_in = Theta3;
          TofHit[1].theta_in = Theta1;
 
          TofHit[0].DT = DT3;
          TofHit[1].DT = DT1;

          TofHit[0].section = 3;
          TofHit[1].section = 1;
        }
        if(TofHit[0].isUorD==2 && TofHit[1].isUorD==1){
          TofHit[0].Esec = E1;
          TofHit[1].Esec = E3;

          TofHit[0].DT = DT1;
          TofHit[1].DT = DT3;

          TofHit[0].theta_in = Theta1;
          TofHit[1].theta_in = Theta3;
 
          TofHit[0].section = 1;
          TofHit[1].section = 3; 
        }
      }

      // *** case 3 *** //
      if(E1!=-1 && E4!=-1){
        if(TofHit[0].isLorR==1 && TofHit[1].isLorR==2){
          TofHit[0].Esec = E1;
          TofHit[1].Esec = E4;

          TofHit[0].theta_in = Theta1;
          TofHit[1].theta_in = Theta4;
          
          TofHit[0].DT = DT1;
          TofHit[1].DT = DT4;

          TofHit[0].section = 1;
          TofHit[1].section = 4;
        }
        if(TofHit[0].isLorR==2 && TofHit[1].isLorR==1){
          TofHit[0].Esec = E4;
          TofHit[1].Esec = E1;

          TofHit[0].theta_in = Theta4;
          TofHit[1].theta_in = Theta1;
           
          TofHit[0].DT = DT4;
          TofHit[1].DT = DT1;

          TofHit[0].section = 4;
          TofHit[1].section = 1;
        }
      }

      // *** case 4 *** //
      if(E2!=-1 && E3!=-1){
        if(TofHit[0].isLorR==1 && TofHit[1].isLorR==2){
          TofHit[0].Esec = E2;
          TofHit[1].Esec = E3;

          TofHit[0].theta_in = Theta2;
          TofHit[1].theta_in = Theta3;
           
          TofHit[0].DT = DT2;
          TofHit[1].DT = DT3;

          TofHit[0].section = 2;
          TofHit[1].section = 3;
        }
        if(TofHit[0].isLorR==2 && TofHit[1].isLorR==1){
          TofHit[0].Esec = E3;
          TofHit[1].Esec = E2;

          TofHit[0].theta_in = Theta3;
          TofHit[1].theta_in = Theta2;
 
          TofHit[0].DT = DT3;
          TofHit[1].DT = DT2;

          TofHit[0].section = 3;
          TofHit[1].section = 2;
        }
      }

      // *** case 5 *** //
      if(E2!=-1 && E4!=-1){
        if(TofHit[0].isUorD==1 && TofHit[1].isUorD==2){
          TofHit[0].Esec = E2;
          TofHit[1].Esec = E4;

          TofHit[0].theta_in = Theta2;
          TofHit[1].theta_in = Theta4;
 
          TofHit[0].DT = DT2;
          TofHit[1].DT = DT4;

          TofHit[0].section = 2;
          TofHit[1].section = 4;
        }
        if(TofHit[0].isUorD==2 && TofHit[1].isUorD==1){
          TofHit[0].Esec = E4;
          TofHit[1].Esec = E2;

          TofHit[0].theta_in = Theta4;
          TofHit[1].theta_in = Theta2;
 
          TofHit[0].DT = DT4;
          TofHit[1].DT = DT2;

          TofHit[0].section = 4;
          TofHit[1].section = 2;
        }
      }

      // *** case 6 *** //
      if(E3!=-1 && E4!=-1){
        if(TofHit[0].isUorD==1 && TofHit[1].isUorD==2){
          TofHit[0].Esec = E3;
          TofHit[1].Esec = E4;

          TofHit[0].theta_in = Theta3;
          TofHit[1].theta_in = Theta4;
 
          TofHit[0].DT = DT3;
          TofHit[1].DT = DT4;

          TofHit[0].section = 3;
          TofHit[1].section = 4;
        }
        if(TofHit[0].isUorD==2 && TofHit[1].isUorD==1){
          TofHit[0].Esec = E4;
          TofHit[1].Esec = E3;

          TofHit[0].theta_in = Theta4;
          TofHit[1].theta_in = Theta3;
 
          TofHit[0].DT = DT4;
          TofHit[1].DT = DT3;

          TofHit[0].section = 4;
          TofHit[1].section = 3;
        }
      }


      // *** spline correction *** //
      for(int i=0; i<2; i++){
        int section = TofHit[i].section;

        double drift_time;
        if(section==1) drift_time = DT1;
        if(section==2) drift_time = DT2;
        if(section==3) drift_time = DT3;
        if(section==4) drift_time = DT4;

        double DT_eval;
        if(section<3) DT_eval=55;
        if(section>2) DT_eval=-55;
        if(section>0){
          TofHit[i].Esec = TofHit[i].Esec / fcorr_z_beta[section-1]->Eval(TofHit[i].beta) * fcorr_z_beta[section-1]->Eval(Beta_norm);

          TofHit[i].Esec = TofHit[i].Esec / fcorr_z_dt[section-1]->Eval(drift_time) * fcorr_z_dt[section-1]->Eval(DT_eval);
        }
      }


      // *** Calculation Theta_out *** //
      double Theta0 = 20.*deg;
      double XA;
      double ZA = 2328.;
      double XC;
      double ZG = 4434.;
      double ZC;
      double XMW3 = -1436.;
      double ZMW3 = 8380;
      double X3lab;
      double Z3lab;
      double Tilt = 14.*deg;
      TVector3 vOut;
      TVector3 vZ = TVector3(0,0,1);
      TVector3 vC;
      TVector3 v3lab;
      for(int i=0; i<2; i++){
        XA = TofHit[i].DT;
        XC = (XA+(ZG-ZA)*tan(TofHit[i].theta_in)) / (1-tan(Tilt)*tan(TofHit[i].theta_in));
        ZC = ZG + XC*tan(Tilt);

        X3lab = TofHit[i].x3*cos(Theta0) + XMW3;
        Z3lab = TofHit[i].x3*sin(Theta0) + ZMW3;

        vC    = TVector3(XC,0,ZC);
        v3lab = TVector3(X3lab,0,Z3lab);
        vOut  = TVector3(X3lab-XC,0,Z3lab-ZC);

        double PathLength = vC.Mag() + vOut.Mag() + 74.;
        PathLength = PathLength/1000.;
        double angle = vZ.Angle(vOut);

        TofHit[i].velocity = PathLength/TofHit[i].tof;
        TofHit[i].beta     = TofHit[i].velocity * m/ns / NPUNITS::c_light;
        TofHit[i].theta_out = angle;
      }

      Z1 = TofHit[0].Esec;
      Z2 = TofHit[1].Esec;


      if(Z1>0 && Z2>0){
        Z1 = fZff_p0 + fZff_p1*Z1 + fZff_p2*Z1*Z1;
        Z2 = fZff_p0 + fZff_p1*Z2 + fZff_p2*Z2*Z2;

        Z1 = sqrt(Z1);
        Z2 = sqrt(Z2);

        Zsum = Z1+Z2;

        iZ1 = (int) round(Z1);
        iZ2 = (int) round(Z2);
        iZsum = iZ1 + iZ2;
      }

      double MagB = 2185*2.2/3584;
      double Leff = 2.067;
      //double rho1 = Leff/abs(2*sin(0.5*(TofHit[0].theta_out - TofHit[0].theta_in)));
      //double rho2 = Leff/abs(2*sin(0.5*(TofHit[1].theta_out - TofHit[1].theta_in)));
      double rho1 = Leff/abs(2*sin(0.5*(TofHit[0].theta_out-TofHit[0].theta_in))*cos(Tilt-0.5*(TofHit[0].theta_out-TofHit[0].theta_in)));
      double rho2 = Leff/abs(2*sin(0.5*(TofHit[1].theta_out-TofHit[1].theta_in))*cos(Tilt-0.5*(TofHit[1].theta_out-TofHit[1].theta_in)));
      double Brho1 = MagB*rho1;
      double Brho2 = MagB*rho2;
      Beta_Z1 = TofHit[0].beta;
      Beta_Z2 = TofHit[1].beta;
      Gamma1 = 1. / sqrt(1 - Beta_Z1 * Beta_Z1);
      Gamma2 = 1. / sqrt(1 - Beta_Z2 * Beta_Z2);

      AoQ1 = Brho1 / (3.10761 * Beta_Z1 * Gamma1);
      AoQ2 = Brho2 / (3.10761 * Beta_Z2 * Gamma2);

      A1 = AoQ1 * iZ1;
      A2 = AoQ2 * iZ2;

      // *** Filling the Fission Fragment Tree *** //
      SofFF->SetTOF(TofHit[0].tof);
      SofFF->SetTOF(TofHit[1].tof);
      SofFF->SetTofPosX(TofHit[0].x);
      SofFF->SetTofPosX(TofHit[1].x);
      SofFF->SetTofPosY(TofHit[0].y);
      SofFF->SetTofPosY(TofHit[1].y);
      SofFF->SetPlastic(TofHit[0].plastic);
      SofFF->SetPlastic(TofHit[1].plastic);
      SofFF->SetPosX3(TofHit[0].x3);
      SofFF->SetPosX3(TofHit[1].x3);
      SofFF->SetThetaIn(TofHit[0].theta_in);
      SofFF->SetThetaIn(TofHit[1].theta_in);
      SofFF->SetThetaOut(TofHit[0].theta_out);
      SofFF->SetThetaOut(TofHit[1].theta_out);


      SofFF->SetBeta(Beta_Z1);
      SofFF->SetBeta(Beta_Z2);
      SofFF->SetGamma(Gamma1);
      SofFF->SetGamma(Gamma2);
      SofFF->SetiZ(iZ1);
      SofFF->SetiZ(iZ2);
      SofFF->SetZ(Z1);
      SofFF->SetZ(Z2);
      SofFF->SetAoQ(AoQ1);
      SofFF->SetAoQ(AoQ2);
      SofFF->SetA(A1);
      SofFF->SetA(A2);
      SofFF->SetBrho(Brho1);
      SofFF->SetBrho(Brho2);

      SofFF->SetDT(TofHit[0].DT);
      SofFF->SetDT(TofHit[1].DT);
      SofFF->SetSection(TofHit[0].section);
      SofFF->SetSection(TofHit[1].section);


      SofFF->SetZsum(Zsum);
      SofFF->SetiZsum(iZsum);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::BeamAnalysis(){
  unsigned int sofsci_size = SofSci->DetectorNbr.size();
  if(sofsci_size==2){
    double beta = SofSci->Beta[0];
    //cout << "Set beta to " << beta << endl;
    SofTrim->SetBeta(beta);
    SofTrim->BuildSimplePhysicalEvent();
    double Zbeam,Theta;
    if(SofTrim->EnergySection.size()>0){
      double Anode1 = SofTrim->EnergySection[0];
      double Anode2 = SofTrim->EnergySection[1];
      double Anode3 = SofTrim->EnergySection[2];
      if(fRunID==13)
        Zbeam = max(Anode1, Anode2);
      else 
        Zbeam = SofTrim->GetMaxEnergySection();

      Theta = SofTrim->Theta[0];

      double TofFromS2    = SofSci->CalTof[0];
      double velocity_mns = SofSci->VelocityMNs[0];
      double Beta         = SofSci->Beta[0];
      double XS2          = SofSci->PosMm[0];
      //double XCC          = SofSci->PosMm[1];
      double XCC=0;
      double YCC=0;
      for(unsigned int i=0; i<SofMwpc->DetectorNbr.size(); i++){
        if(SofMwpc->DetectorNbr[i]==1){
          XCC = SofMwpc->PositionX1[i];
          YCC = SofMwpc->PositionY[i];
        }
      }

      double LS2;
      LS2 = fLS2_0;//*(1 + fK_LS2*Theta);
      velocity_mns = LS2/TofFromS2;
      Beta = velocity_mns * m/ns / NPUNITS::c_light;
      double Gamma        = 1./(TMath::Sqrt(1 - TMath::Power(Beta,2)));
      double Brho = fBrho0 * (1 - XS2/fDS2 - XCC/fDCC);
      double AoQ  = Brho / (3.10716*Gamma*Beta);

      // Y dependence correction //
      double Y_p0 = 23943.8;
      double Y_p1 = 12.362;
      Zbeam = Zbeam/(Y_p0 + Y_p1*YCC)*Y_p0;


      // Z calibration //
      Zbeam = fZbeam_p0 + fZbeam_p1*Zbeam + fZbeam_p2*Zbeam*Zbeam;
      Zbeam = sqrt(Zbeam);
      double A = AoQ * round(Zbeam);

      // Last beta correction //
      double Beta_norm = 0.8355;
      Zbeam = Zbeam/(fZBeta_p0 + fZBeta_p1*Beta)*(fZBeta_p0 + fZBeta_p1*Beta_norm);

      // Filling Beam tree
      SofBeamID->SetZbeam(Zbeam);
      SofBeamID->SetAoQ(AoQ);
      SofBeamID->SetAbeam(A);
      SofBeamID->SetBeta(Beta);
      SofBeamID->SetGamma(Gamma);
      SofBeamID->SetBrho(Brho);
      SofBeamID->SetXS2(XS2);
      SofBeamID->SetXCC(XCC);
      SofBeamID->SetYCC(YCC);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
void Analysis::LoadSpline(){
  TString input_path = "./calibration/SofTwim/spline/";

  TString rootfile = input_path + "spline_beta.root";
  TFile* ifile = new TFile(rootfile,"read");

  TString splinename;
  if(ifile->IsOpen()){
    cout << "Loading Beta spline for fission fragment analysis..." << endl;
    for(int i=0; i<4; i++){
      splinename = Form("spline_beta_sec%i",i+1);
      fcorr_z_beta[i] = (TSpline3*) ifile->FindObjectAny(splinename);
    }
    ifile->Close();
  }
  else
    cout << "File " << rootfile << " not found!" << endl;

  //*** ***//
  rootfile = input_path + "spline_dt.root";
  ifile = new TFile(rootfile,"read");

  if(ifile->IsOpen()){
    cout << "Loading DT spline for fission fragment analysis..." << endl;
    for(int i=0; i<4; i++){
      splinename = Form("spline_dt_sec%i",i+1);
      fcorr_z_dt[i] = (TSpline3*) ifile->FindObjectAny(splinename);
    }
    ifile->Close();
  }
  else
    cout << "File " << rootfile << " not found!" << endl;

}


////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitParameter(){
  fLS2_0 = 135.614;
  fDS2   = 8000;
  fDCC   = -10000;
  fK_LS2 = -30e-8;

  fBrho0 = 12.3255;
  fRunID = 5;

  // Beam parameter //
  fZBeta_p0 = 1;
  fZBeta_p1 = 0;
  fZbeam_p0 = 1651.57;
  fZbeam_p1 = 0.0876127;
  fZbeam_p2 = 4.02563e-6;

  // FF parameter //
  fZff_p0 = 2.80063;
  fZff_p1 = 6.91985e-2;
  fZff_p2 = 1.01598e-7;

  if(fRunID==1 || fRunID==2){
    //fBrho0 = 10.6813; // 180Hg
    fBrho0 = 10.6955; // 180Hg
    fZbeam_p0 = -5303.06;
    fZbeam_p1 = 0.674945;
    fZbeam_p2 = -8.32085e-6;

    fZBeta_p0 = 72.946;
    fZBeta_p1 = 6.0644;
  }
  if(fRunID==3){
    fBrho0 = 10.8183; // 182Hg
    fZbeam_p0 = -2737.25;
    fZbeam_p1 = 0.452017;
    fZbeam_p2 = -3.48831e-6;

    fZBeta_p0 = 76.6738;
    fZBeta_p1 = 1.60128;
  }
  if(fRunID==4){
    fBrho0 = 10.9558; // 184Hg
    fZbeam_p0 = -5044.61;
    fZbeam_p1 = 0.639986;
    fZbeam_p2 = -7.3077e-6;
  }
  if(fRunID==5){
    fBrho0 = 10.8138; // 187Pb
    fZbeam_p0 = -2858.72;
    fZbeam_p1 = 0.454064;
    fZbeam_p2 = -3.36443e-6;

    fZBeta_p0 = 71.0975;
    fZBeta_p1 = 10.7007;
  }
  if(fRunID==6){
    fBrho0 = 10.9476; // 189Pb
    fZbeam_p0 = 1590.66;
    fZbeam_p1 = 0.0956455;
    fZbeam_p2 = 3.84585e-6;

    fZBeta_p0 = 74.6063;
    fZBeta_p1 = 6.4635;
  }
  if(fRunID==7){
    fBrho0 = 10.6814; // 175Pt
    fZbeam_p0 = 459.68;
    fZbeam_p1 = 0.162277;
    fZbeam_p2 = 3.10164e-6;

    fZBeta_p0 = 66.9433;
    fZBeta_p1 = 10.8664;
  }
  if(fRunID==8){
    fBrho0 = 11.0864; // 204Fr
    fZbeam_p0 = 4122.94;
    fZbeam_p1 = -0.119867;
    fZbeam_p2 = 8.29115e-6;

    fZBeta_p0 = 63.9575;
    fZBeta_p1 = 25.1988;
  }
  if(fRunID==9){
    fBrho0 = 11.2712; // 207Fr
    fZbeam_p0 = -1752.27;
    fZbeam_p1 = 0.346018;
    fZbeam_p2 = -8.64673e-7;

    fZBeta_p0 = 63.9575;
    fZBeta_p1 = 25.1988;
  }
  if(fRunID==10){
    fBrho0 = 11.0955; // 199At run 423 & 424
    fZbeam_p0 = -116.425;
    fZbeam_p1 = 0.218256;
    fZbeam_p2 = 1.62399e-6;

    fZBeta_p0 = 61.3889;
    fZBeta_p1 = 25.8908;
  }
  if(fRunID==11){
    fBrho0 = 10.9970; // 199At run 425 & 426
    fZbeam_p0 = -116.425;
    fZbeam_p1 = 0.218256;
    fZbeam_p2 = 1.62399e-6;

    fZBeta_p0 = 61.3889;
    fZBeta_p1 = 25.8908;
  }
  if(fRunID==12){
    fBrho0 = 10.8697; //197At
    fZbeam_p0 = -2683.52;
    fZbeam_p1 = 0.422551;
    fZbeam_p2 = -2.44471e-6;

    fZBeta_p0 = 62.9188;
    fZBeta_p1 = 22.8827;
  }
  if(fRunID==13){
    fBrho0 = 11.3418; // 216Th
    fZbeam_p0 = 1651.57;
    fZbeam_p1 = 0.0876127;
    fZbeam_p2 = 4.02563e-6;

    fZBeta_p0 = 38.879;
    fZBeta_p1 = 58.7667;
  }
  if(fRunID==14){
    fBrho0 = 11.5067; // 221Pa
    fZbeam_p0 = 186.892;
    fZbeam_p1 = 0.20739;
    fZbeam_p2 = 1.61797e-6;
  }
  if(fRunID==15){
    fBrho0 = 12.3352;
  }

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  SofBeamID->Clear();
  SofFF->Clear();
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  RootOutput::getInstance()->GetTree()->Branch("RunID",&RunID,"RunID/I");
  RootOutput::getInstance()->GetTree()->Branch("SofBeamID","TSofBeamID",&SofBeamID);
  RootOutput::getInstance()->GetTree()->Branch("SofFissionFragment","TSofFissionFragment",&SofFF);

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


