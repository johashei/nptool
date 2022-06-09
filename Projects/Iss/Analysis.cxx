/*****************************************************************************
 * Copyright (C) 2009-2019    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: M Labiche  contact address: marc.labiche@stfc.ac.uk      *
 *                                                                           *
 * Creation Date  : july 2019                                                *
 * Last update    :                                                          *
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
using namespace std;
#include"Analysis.h"
#include"NPAnalysisFactory.h"
#include"NPDetectorManager.h"
#include"NPOptionManager.h"
#include"NPFunction.h"
#include"RootOutput.h"
#include"RootInput.h"
#include"NPPhysicalConstants.h"


////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  Iss = (TIssPhysics*)  m_DetectorManager->GetDetector("ISSDet");

    
  myInit= new TInitialConditions();
  myIntCoord= new TInteractionCoordinates();
  
  InitOutputBranch();
  InitInputBranch();


  Bfield= Iss->GetNominalMField() ;
  cout << "Nominal Magnetic field= " <<  Bfield  << " milli(?) Tesla" << endl;
  Bfield=Bfield*1000.;
  cout << "Nominal Magnetic field= " <<  Bfield  << " Tesla" << endl;

  //LightCD2 = EnergyLoss("Deuteron_CD2.G4table","G4Table",100 );
  //LightCD2 = EnergyLoss("He3_CD2.G4table","G4Table",100 );
  LightCD2 = EnergyLoss("proton_CD2.G4table","G4Table",100);
  //LightAl = EnergyLoss("proton_Al.G4table","G4Table",100);
  //LightSi = EnergyLoss("proton_Si.G4table","G4Table",100);
  //BeamCD2 = EnergyLoss("Sn132_CD2.G4table","G4Table",100);
  BeamCD2 = EnergyLoss("Mg28_CD2.G4table","G4Table",100);
  //BeamCD2 = EnergyLoss("Ce146_CD2.G4table","G4Table",100);
  //BeamCD2 = EnergyLoss("Rn212_CD2.G4table","G4Table",100);
  //BeamCD2 = EnergyLoss("Pb190_CD2.G4table","G4Table",100);
  //BeamCD2 = EnergyLoss("Sn108_CD2.G4table","G4Table",100);
  //BeamCD2 = EnergyLoss("Ne17_CD2.G4table","G4Table",100);
  //BeamCD2 = EnergyLoss("Hg206_CD2.G4table","G4Table",100);

  myReaction = new NPL::Reaction();
  myReaction->ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());

   TargetThickness = m_DetectorManager->GetTargetThickness();  // * by micrometer to convert in mm  
   cout << "Target Thickness = " <<  TargetThickness << " mm" << endl; 


   OriginalBeamEnergy = myReaction->GetBeamEnergy();
   cout << "Beam Energy before target= " <<  OriginalBeamEnergy << " MeV" << endl; 
 

  pCharge= myReaction->GetNucleus3()->GetZ();

  //double pAtomWeight= myReaction->GetNucleus3()->GetA();
  //double pMassExcess= myReaction->GetNucleus3()->GetMassExcess();
  //pMass= (pAtomWeight*amu_c2 + pMassExcess/1000)/amu_c2;
  //cout << "Mass of light particle= " <<  pMass << endl; 
  //pMass= (myReaction->GetNucleus3().Mass())/amu_c2;
  //pMass=1.00782503207;  // MeV/c2 - for 132Sn(d,p)
  //pMass=2.0141017778;  // MeV/c2 - for 224Ra(d,d)
  //pMass=12*931.49;  // MeV/c2 - for 224Ra(C,C)
  pMass1= (myReaction->GetNucleus1()->Mass())/amu_c2; // HI beam
  pMass2= (myReaction->GetNucleus2()->Mass())/amu_c2; // tgt
  pMass3= (myReaction->GetNucleus3()->Mass())/amu_c2; // light ejectile
  pMass4= (myReaction->GetNucleus4()->Mass())/amu_c2; // beam-like


  cout << "Charge of light particle= " <<  pCharge << endl;
  cout << "Atomic Mass unit= " <<  amu_c2  << endl; 
  cout << "Mass of light particle= " <<   myReaction->GetNucleus3()->Mass() << endl; 
  cout << "Mass of light particle= " <<  pMass3 << " amu" << endl;
  cout << "Mass of  proj= " <<  pMass1 << " amu" << endl; 
  cout << "Mass of  targ= " <<  pMass2 << " amu" << endl; 
  cout << "Mass of heavy particle= " <<  pMass4 << " amu" << endl;  
  cout << "Reaction QValue= " <<  myReaction->GetQValue() << " MeV" <<endl; 


   Rand = new TRandom3();
   DetectorNumber = 0 ;
   ThetaNormalTarget = 0 ;

   ThetaIssSurface = 0;
   X_Iss = 0 ;
   Y_Iss = 0 ;
   Z_Iss = 0 ;
   Si_E_Iss = 0 ;
   E_Iss = 0;
   Si_X_Iss = 0;
   Si_Y_Iss = 0;

  double BeamEnergy = BeamCD2.Slow(OriginalBeamEnergy,TargetThickness*0.5,0);
  myReaction->SetBeamEnergy(BeamEnergy);
  cout << "Beam energy set at " << BeamEnergy << " MeV" << endl;  
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  // Reinitiate calculated variable
  ReInitValue();

  double XTarget = -1001;
  double YTarget = -1001;

  TVector3 BeamDirection = TVector3(0,0,1);

  double BeamEnergy = BeamCD2.Slow(OriginalBeamEnergy,TargetThickness/2,0);
  myReaction->SetBeamEnergy(BeamEnergy); // set the new beam energy taking into account the beam energy loss in target befor reaction
 
   double Theta_extrap=0. ;
   double ThetaDet =0. ; 
   //double ExcitationE=0. ; double ExcitationE_extrap=0. ; 


   TVector3 Pos; // for Iss


   XTarget=myInit->GetIncidentPositionX();
   YTarget=myInit->GetIncidentPositionY();

 
  //if(Iss->Strip_E.size()>0)cout << "multiplicity= " << Iss->Strip_E.size()<< endl;


  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////// LOOP on ISS //////////////////
  if(Iss->Strip_E.size()>=1){

 
    /************************************************/
    // Part 1 : Interaction position in detector

    if(XTarget>-1000 && YTarget>-1000){ // XTarget and YTarget can be known experimentaly only if there are beam tracker detectors !!!

      //Pos = Iss->GetPositionOfInteraction(0,1); // Lab Position of the light particle interaction in the Si Detector for Iss (random inside the pixel)


      Pos = Iss->GetPositionOfInteraction(0,0); // Lab Position of the light particle interaction in the Si Detector for Iss  (center of pixel)


	  X_meas = Pos(0)/mm;
	  Y_meas = Pos(1)/mm;
	  Z_meas = Pos(2)/mm;

      //Z_Iss = myIntCoord->GetDetectedPositionZ(0);

    }
    
    else{
      BeamDirection = TVector3(-1000,-1000,-1000);
    }
    
    /************************************************/
    
    /************************************************/
    // Part 2 : Impact Energy
    Energy = ISS_EDep = 0;
    //Energy = Iss->GetEnergyDeposit();
    if(Iss->Strip_E[0]>0){
         Energy = Iss->Strip_E[0];
    }

    ISS_DetId= Iss->DetectorNumber[0];  // Assuming Multiplicity 1
      
	// Detector intrinsic energy resolution
    //ISS_EDep   = Energy; // resolution already taken into account in simulation output for Is !!
    //ISS_EDep   = Rand->Gaus(Energy, 0.070/2.35); // 70 keV sigma resolution for Iss !!
    ISS_EDep   = Rand->Gaus(Energy, 0.060/2.35); // 60 keV sigma resolution for Iss !!
    //ISS_EDep   = Rand->Gaus(Energy, 0.020/2.35); // 20 keV sigma resolution for Iss !!
    //ISS_EDep   = Rand->Gaus(Energy, 0.040/2.35); // 40 keV sigma resolution for Iss !!
    //ISS_EDep   = Rand->Gaus(Energy, 0.080/2.35); // 80 keV sigma resolution for Iss !!
    //ISS_EDep   = Rand->Gaus(Energy, 0.10/2.35); // 100 keV sigma resolution for Iss !!

    ISS_Time= Iss->Strip_T[0]; 

    /************************************************/
    
    /************************************************/

    // Part 3 for Iss: Z and theta extrapolation
    ZonBeamAxis(Pos,ISS_EDep,pCharge,pMass3,Bfield,&Z_extrap,&ThetaDet,&Theta_extrap,&Nb_Iter);
    if(Theta_extrap>=0){ // Todo: check in ZonBeamAxis why sometime nothing is returned for Theta_extrap !! 
      //cout << "Theta_extrap:" << Theta_extrap << endl;

     ThetaLab=Theta_extrap; // in degrees

   // Target Correction (ELab= EDep + Energy loss in target given ThetaLab )
      ISS_ELab = ISS_EDep; // without target correction     
      //cout << "ELab1=" << ELab <<endl;
      //cout << "ThetaLab:" << ThetaLab << endl;

      ISS_ELab = LightCD2.EvaluateInitialEnergy( ISS_ELab ,TargetThickness/2., ThetaLab*deg); // with target correction but neglecting angular straggling in target


    /************************************************/
    
    /************************************************/
    // Part 4 : Excitation Energy Calculation
      Ex = myReaction -> ReconstructRelativistic( ISS_ELab , ThetaLab*deg );
   
    /************************************************/
    
    /************************************************/
    // Part 5 : Theta CM Calculation
      ThetaCM  = myReaction -> EnergyLabToThetaCM( ISS_ELab , ThetaLab*deg )/deg;    
    
    /************************************************/
	// Part 6 : Qvalue
	
   // !!!! Reaction dependent parameters for QValue !!!! :

   // initialisation:
     //parA=1.0076 ; // (MeV) For132Sn(d,p) at 10 MeV/u (independent of B)
     //parB=-9.9247 ; // For132Sn(d,p) at 10 MeV/u (independant on B)
     //parC=-.1041 ; // (cm) For132Sn(d,p) at 10 MeV/u and B=1.5 T

     //parC=-.139 ; // (cm) For132Sn(d,p) at 10 MeV/u and B=2 T
     //parC=-.2086 ; // (cm) For132Sn(d,p) at 10 MeV/u and B=3 T

     //parA=1.0089 ; // (MeV) For 224Ra(d,d) at 4 MeV/u
     //parB=0.0001 ; // For 224Ra(d,d) at 4 MeV/u
     //parC=-.1061 ; // (cm) For 224Ra(d,p) at 4 MeV/u
     // parA= 1.0536 ; // (MeV) For 224Ra(C,C) at 4.5 MeV/u
     // parB= 0.0008 ; // For 224Ra(C,C) at 4.5 MeV/u
     // parC= -0.5628 ; // (cm) For 224Ra(C,C) at 4.5 MeV/u

    /* Calculating parA, parb, and parC : */
	 double Vcm;
	 double Vbeam;
	 double gbeam;
	 //cout << "%%%%%% " << endl;
			 
		Vcm= sqrt((2*myReaction->GetBeamEnergy()*(pMass1*amu_c2))) / ((pMass1+pMass2)*amu_c2) ; // in unit of c  & non relativistic
		Vbeam = sqrt(2*myReaction->GetBeamEnergy()/(pMass1*amu_c2)); // in unit of c - non relativist
		//cout << "Vbeam=" <<Vbeam << endl;
		//cout << "Vcm=" <<Vcm << endl;	    

  	    //cout << Vcm  << " " << Vbeam << endl;
		parA= (pMass3+pMass4)/pMass4;
		//parB= (parA*((pMass3/2.)*amu_c2)*pow(Vcm,2)  -  ((pMass1+pMass2)/pMass2)*((pMass1/2.)*amu_c2)*(pow(Vbeam-Vcm,2)));
		parB= parA*((pMass3/2.)*amu_c2)*pow(Vcm,2)  -  myReaction->GetBeamEnergy()*(pMass2/(pMass1+pMass2));
		parC=-(parA * pMass3 *amu_c2* Vcm) / ((c_light/(cm/s)) *2*pi*pMass3*1.660539e-27/(Bfield*pCharge*1.60218e-19)) ;
		//parC=-(parA * pMass3 *amu_c2* Vcm) / ((c_light/(cm/s)) *54.6667e-9) ;

    /* Calculating Q value non relativistic*/
	  QValue= parA*ISS_ELab + parB + parC*(Z_extrap/cm);  // Peter's parameters (Z_extrap/cm to transform Z_extrap from mm to cm)
	//if(Iss->GetTimeDetected()<31)cout << "Qv= " <<  parA*ELab + parB + parC*(Z_extrap/cm) << endl;
	  //QValue= QValue + 6.44808e-05*(Z_extrap)-0.00985694;  // Peter's parameters (Z_extrap/cm to transform Z_extrap from mm to cm)
	//cout << "QValue=" << QValue << endl;

    // or using relativistic reconstruction of Ex:
	  QValueR= ((pMass1+pMass2)-(pMass3+pMass4))*amu_c2 + Ex;

    /* Calculating maximum distance to beam axis*/
   	  MaxDAxis(ISS_ELab,pCharge,pMass3,Bfield,&Daxis);
   }

  }//end loop Iss



}


////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
			cout << "Ebeam half way through Target= " << BeamEnergy << endl;
			cout << "Tcycl= " << 65.6*( int(pMass3+0.5)/(pCharge*Bfield)) << "ns" << endl;
			cout << "Tcycl= " << 2*pi*pMass3*1.660539e-27/(Bfield*pCharge*1.60218e-19) << "s" << endl;
			cout << "pMass1= " <<pMass1 << " pMass4=" << pMass4  << endl;
			cout << "pMass2= " <<pMass2 << " pMass3=" << pMass3  << endl;
			cout << "parA= " <<parA << endl;
			cout << "parB= " <<parB << endl;
			cout << "parC= " <<parC << endl;
	
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("ISS_EDep",&ISS_EDep,"ISS_EDep/D");
  RootOutput::getInstance()->GetTree()->Branch("ISS_ELab",&ISS_ELab,"ISS_ELab/D");
  RootOutput::getInstance()->GetTree()->Branch("ISS_DetId",&ISS_DetId,"ISS_DetId/I");
  RootOutput::getInstance()->GetTree()->Branch("ISS_Time",&ISS_Time,"ISS_Time/D");
  RootOutput::getInstance()->GetTree()->Branch("Z_Iss",&Z_extrap,"Z_Iss/D"); // projected to beam axis
//  RootOutput::getInstance()->GetTree()->Branch("Zextrap",&Z_extrap,"Z_extrap/D"); // projected to beam axis

  RootOutput::getInstance()->GetTree()->Branch("Ex",&Ex,"Ex/D");
  RootOutput::getInstance()->GetTree()->Branch("Qv",&QValue,"QValue/D");
  RootOutput::getInstance()->GetTree()->Branch("QvR",&QValueR,"QValueR/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab,"ThetaLab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM",&ThetaCM,"ThetaCM/D");
  RootOutput::getInstance()->GetTree()->Branch("Xmeas",&X_meas,"X_meas/D");  // measured
  RootOutput::getInstance()->GetTree()->Branch("Ymeas",&Y_meas,"Y_meas/D");  // measured
  RootOutput::getInstance()->GetTree()->Branch("Zmeas",&Z_meas,"Z_meas/D");  // measured
  RootOutput::getInstance()->GetTree()->Branch("NbIter",&Nb_Iter,"Nb_Iter/I"); // 
  RootOutput::getInstance()->GetTree()->Branch("Daxis",&Daxis,"Daxis/D"); // 

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  //RootInput:: getInstance()->GetChain()->SetBranchAddress("InitialConditions",&myInit );
  //RootInput:: getInstance()->GetChain()->SetBranchAddress("InteractionCoordinates",&myIntCoord );
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  ISS_ELab = -100;
  ISS_EDep = -100;
  ISS_DetId = -10;
  ISS_Time = -10;
  Z_Iss= -1000; 
  Ex= -100;
  QValue = -10 ;
  QValue = -10 ;
  ThetaLab = -100;
  ThetaCM = -100;
  X_meas = -1000;
  Y_meas = -1000;
  Z_meas = -1000;
  Z_extrap = -1000;
  Nb_Iter = -10;
  Daxis = -0.0001;

 
  E_S1=S1ELab=-100;
  Time_S1=-100;
  Sector_S1=-10;
  Ring_S1=-10;

  //S1E=-100;
  S1Ring=-10;
  S1Sector=-10;
  S1Time=-10;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ZonBeamAxis (TVector3 A, double nrj , int q, double Mp, double B, double *Z_extrap, double *ThetaDet, double *Theta_extrap, int *Nb_Iter)
{
  //   double Z = acos( (A.Dot(B)) / (A.Mag()*B.Mag()) );
  double Zinit,Zf1, Zf2;
  double Zeta_est=0;  // = Theta lab of emitted particle
  double Rho_max=0; 
  double Rprim, alpha;
  double DeltaZ1,DeltaZ2;
  int Nb_iteration;
  //cout << "############################################################"<< endl;
  //cout << "nrj=" << nrj << " q=" << q << " Mp=" << Mp << " B=" << B << endl;

  Zf1=Zinit= A(2)/m; // change Zdet from mm to m 
  //cout << "Zf1=" << Zf1 << endl;

  Rprim=sqrt(pow(A(0),2)+pow(A(1),2))/m; // from mm to m 
  //cout << "x=" << A(0) << "mm" << " y=" << A(1) << "mm"<< endl;
  //cout << "Rprim=" << Rprim << "m"<< endl;

  Zeta_est= acos((Zf1*q*B)/(0.9046*sqrt(nrj*Mp)));  
  //cout << "cos_Zeta_est=" << (Zf1*q*B)/(0.9046*sqrt(nrj*Mp)) << " Zeta_est=" << Zeta_est*180/3.14159 << "deg"<< endl;

  *ThetaDet=Zeta_est/deg; // first estimation of Theta

  //cout << "ThetaDet=" <<  Zeta_est/deg  << "deg"<< endl;

  Rho_max=0.2879*sqrt(nrj*Mp)*(sin(Zeta_est))/(q*B);
  //cout << " Rho_max=" << Rho_max << endl;

    alpha=2*asin(Rprim/Rho_max);
    //cout << "  alpha=" << alpha*180/3.14152 << "deg" << endl;

     DeltaZ1=Zf1*alpha/(twopi);
     //cout << "   DeltaZ1= " << DeltaZ1 << endl;

      Zf2=DeltaZ1+Zf1;
      //cout << "   Zf2= " << Zf2 << endl;

	Nb_iteration=1;

// Let's iterate:

       DeltaZ2=DeltaZ1;  // initialisation of DeltaZ2

       //while( sqrt((Zf2-Zf1)*(Zf2-Zf1)) > 0.000005 && DeltaZ1*DeltaZ2>0 )  // we itirate until we get a precision better than 0.005 mm OR until DeltaZ change sign (<=> to alpha changing sign).
     while( sqrt(pow((Zf2-Zf1),2))*m > 0.00005 )  // we itirate until we get a precision better than 0.00005 mm 
    {
       Zf1=Zf2; 

        Zeta_est= acos((Zf1*q*B)/(0.9046*sqrt(nrj*Mp)));  // estimating the new theta lab
        //cout << "cos_Zeta_est=" << (Zf1*q*B)/(0.9046*sqrt(nrj*Mp)) << " theta_est=" << Zeta_est*180/3.14159 << "deg" << endl;

         Rho_max=0.2879*sqrt(nrj*Mp)*(sin(Zeta_est))/(q*B);  // estimating the new Rho_max
         //cout << " Rho_max=" << Rho_max << endl;
 
         //alpha=alpha-2*asin(Rprim/Rho_max);  // estimating the new alpha
         alpha=2*asin(Rprim/Rho_max); // estimating the new alpha with the new Rho_max
          //cout << "  alpha=" << alpha*180/3.14159 << " deg" << endl;

          //DeltaZ2=Zf1*alpha/(2*3.14159);
          DeltaZ2=Zinit*alpha/(twopi); 
          //cout << "   DeltaZ2= " << DeltaZ2 << endl;

         //Zf2=DeltaZ2 + Zf1;
          Zf2=DeltaZ2 + Zinit;  // estimating the new Zf2
         //cout << "   New Zf2= " << Zf2 << endl;

     
		Nb_iteration++;
    }
    
    //cout << "Nb_iteration= " << Nb_iteration<< endl;

      *Z_extrap=Zf2*m;   // retuns in mm
      *Theta_extrap=Zeta_est/deg;  // returns Theta_extrap in deg
	  *Nb_Iter=Nb_iteration;
      return; // returning a value in mm
}


////////////////////////////////////////////////////////////////////////////////
// To calculate the maximum distance to axis
void Analysis::MaxDAxis (double nrj , int q, double Mp, double B, double *D)
{
	double Vt; // transversal velocity
	double Rc; //  rayon de courbure
	double Result; //
	
	
	//cout << " q= " << q << endl;
	//cout << " Mp= " << Mp << endl;
	
	Vt= sqrt(2*nrj/(Mp*amu_c2))*c_light; // in m/s
	
	Rc= Mp*1.660539e-27*Vt/(q*1.60218e-19*B); // in m
	
	//cout << "nrj=" << nrj << endl;
	//cout << "Vt=" << Vt << endl;
	//cout << "q=" << q << endl;
	//cout << "B=" << B << endl;
	//cout << "Rc=" << Rc << endl;
	
	Result= 2.*Rc*100.; // * by 100. to convert in cm
	
	//cout << "Result=" << Result << endl;
	
	*D = Result;  
	return;
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
class proxy{
  public:
    proxy(){
      NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct);
    }
};

proxy p;
}

