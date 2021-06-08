#include "MinimizeBeamSpot.h"

double devE(const double* parameter){
  //Beam energy: 7.7 [MeV/A] * 47 [A] = 361.9 [MeV]
  static NPL::Reaction reaction("47K(d,p)48K@362");

  //Beam spot offset
  TVector3 offset(parameter[0],parameter[1],parameter[2]);

  //Other variable initilizations
  unsigned int size = pos.size();
  TVector3 MugastNormal;
  double dE,Theta;
  TVector3 dir;

  //Initilize histogram
  //canv->Clear();
  //canv->ResetDrawn();
  h->Reset();
  h1->Reset();
  h2->Reset();
  h3->Reset();
  h4->Reset();
  h5->Reset();
  h7->Reset();

  double FitResultMatrix[7][5];

  //Loop over events
  for(unsigned int i = 0 ; i < size ; i++){
    //Particle path vector
    dir=*(pos[i])-offset;

    //Define normal vector for the MG# of detection
    DetectorSwitch(detnum[i], MugastNormal);

    //Angle leaving target, angle entering MUGAST & energy deposited in MUGAST
    double ThetaTarget = dir.Angle(TVector3(0,0,1));
    double ThetaMugast = dir.Angle(MugastNormal);
    double Energy = energy[i];

    //Energy loss in Al
    Energy=Al.EvaluateInitialEnergy(
		Energy,                      //energy after Al 
		0.4*micrometer,              //thickness of Al
		ThetaMugast);                //angle impinging on MUGAST
    //Energy loss in target
    Energy=CD2.EvaluateInitialEnergy(
		Energy,                      //energy after leaving target
		0.5*parameter[3]*micrometer, //pass through half target
		ThetaTarget);                //angle leaving target

    //Final value of Ex
    double Ex = reaction.ReconstructRelativistic(Energy,ThetaTarget);
    
    //Fill histograms with Ex
    h->Fill(Ex);
    DetectorSwitch(detnum[i], Ex);
  }
  
  //Initilise, Draw & Fit histograms
  InitiliseCanvas(FitResultMatrix);
  
  //Write vals to screen
  if(flagDraw){cout << "==================================================" << endl;}
  cout << "Mean: "     << FitResultMatrix[mgSelect][0]
	               << " +/- "
		       << FitResultMatrix[mgSelect][1]
       << "    StdDev: " << FitResultMatrix[mgSelect][2] 
	               << " +/- "
		       << FitResultMatrix[mgSelect][3]
       << "    Thick: " << parameter[3] << " um"
       << "    Fit Chi2/NDF = " << FitResultMatrix[mgSelect][4]
       << endl;

  //Adapt the metric as needed
  return sqrt( pow(FitResultMatrix[mgSelect][0]-refE,2) + pow(0.1*FitResultMatrix[mgSelect][2],2) );
}
////////////////////////////////////////////////////////////////////////////////
void MinimizeBeamSpot(){

  //Read data in
  LoadFile();

  //Output formatting
  cout << fixed << showpoint << setprecision(6) << showpos;
  
  //Read in
  cout << "==================================================" << endl;
  cout << "=--------- SELECT TELESCOPE TO MINIMIZE ---------=" << endl;
  cout << "= Type MG# of telescope metric to use, or type 0 =" << endl;
  cout << "= to use the sum of all MG's                     =" << endl;
  cout << "==================================================" << endl;
      cin >> mgSelect;
      if(mgSelect==7){mgSelect=6;} // Correct the input for MG7
  cout << "==================================================" << endl;

  //Start with beam (0,0,0) and 4.76um 0.5mg/c2 target
  double parameter[4] = {0.0, 0.0, 0.0, 4.76};   
  devE(parameter);

  //Function with 4 parameter XYZ and Target thickness
  auto func = ROOT::Math::Functor(&devE,4);
 
  //Initilise minimizer
  auto minim = ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"); 
  minim->SetPrintLevel(0);
  minim->SetPrecision(1e-10); 

  //Set minimizer function
  minim->SetFunction(func);

  //Assign variable limits
  minim->SetLimitedVariable(0,"X",parameter[0],0.01,-10,10);
  minim->SetLimitedVariable(1,"Y",parameter[1],0.01,-10,10);
  minim->SetLimitedVariable(2,"Z",parameter[2],0.01,-5,5);
  minim->SetLimitedVariable(3,"T",parameter[3],0.01,4.76*0.5,4.76*1.5);

  //Don't draw iterations of minimizer
  flagDraw = 0;
  //canv->SetBatch(kTRUE);
  gROOT->SetBatch(kTRUE);
  
  //Shrink it, babeyyy
  minim->Minimize(); 
  
  //Draw minimal value
  flagDraw = 1;
  gROOT->SetBatch(kFALSE);

  //Pull values from minimizer
  const double* x = minim->X();
  cout << "==================================================" << endl;
  cout << "=---------------- FINAL PEAK FITS ---------------=" << endl;
  cout << "==================================================" << endl;
    devE(x);
//    canv->DrawClone();
  cout << "==================================================" << endl;
  cout << "=------------ RESULTS OF MINIMIZATION -----------=" << endl;
    if(mgSelect==6){mgSelect=7;} // Correct the input for MG7
  cout << "=------------------- USING MG " << mgSelect << " -----------------=" << endl;
  cout << "==================================================" << endl;
  cout << "\t\tX = " << x[0] << " mm" << endl;
  cout << "\t\tY = " << x[1] << " mm" << endl;
  cout << "\t\tZ = " << x[2] << " mm" << endl;
  cout << "\t\tT = " << x[3] << " um" << endl;
  cout << "==================================================" << endl;
}
