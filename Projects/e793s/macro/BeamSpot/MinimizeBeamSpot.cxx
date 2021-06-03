#include "MinimizeBeamSpot.h"

double devE(const double* parameter){
  //Beam energy: 7.7 [MeV/A] * 47 [A] = 361.9 [MeV]
  static NPL::Reaction reaction("47K(d,p)48K@362");

  //Beam spot offset
  TVector3 offset(parameter[0],parameter[1],parameter[2]);
  unsigned int size = pos.size();

  double dE,Theta;
  TVector3 dir;

  //Initilize histogram
  h->Reset();
  h1->Reset();
  h2->Reset();
  h3->Reset();
  h4->Reset();
  h5->Reset();
  h7->Reset();

  //Loop over events
  for(unsigned int i = 0 ; i < size ; i++){
    //Particle path vector
    dir=*(pos[i])-offset;

    //Detected energy, and angle of particle leaving target
    double Theta= dir.Angle(TVector3(0,0,1));
    double Energy = energy[i];

    //NOTE!!! Not calucualting energy loss in Al???
    //Energy loss in target
    Energy=CD2.EvaluateInitialEnergy(
		    Energy,                      //energy after leaving target
		    0.5*parameter[4]*micrometer, //pass through half target
		    Theta);                      //angle leaving target

    //Final value of Ex
    double Ex = reaction.ReconstructRelativistic(Energy,Theta);
    
    //Fill histogram with ERROR in Ex!
    //h->Fill(Ex-refE); 
    h->Fill(Ex); 

    switch(detnum[i]){
      case 1:
        h1->Fill(Ex); 
        break;
      case 2:
        h2->Fill(Ex); 
        break;
      case 3:
        h3->Fill(Ex); 
        break;
      case 4:
        h4->Fill(Ex); 
        break;
      case 5:
        h5->Fill(Ex); 
        break;
      case 7:
        h7->Fill(Ex); 
        break;
      default:
        cout << "ERROR! Invalid detnum: " << detnum[i] << " @" << i << endl;
        return 1; // Exit code
    }

  }
  //End loop over events

  //Write vals to screen
  cout << "Mean: " << h->GetMean() 
       << "\t StdDev: " << h->GetStdDev() 
       << "\t Thickness??: " << parameter[4] 
       << endl;

  //Draw histogram(s)
  h->Draw();
  if(flagDraw){ InitiliseCanvas(); }

  //Adapt the metric as needed
  return sqrt( pow(h->GetMean()-refE,2) + pow(0.1*h->GetStdDev(),2) );
}
////////////////////////////////////////////////////////////////////////////////
void MinimizeBeamSpot(){

  // Read data in
  LoadFile();

  // Start with beam (0,0,0) and 4.7um 0.5mg/c2 target
  double parameter[4] = {0.0, 0.0, 0.0, 4.7};   
  devE(parameter);

  // Function with 4 parameter XYZ and Target thickness
  auto func = ROOT::Math::Functor(&devE,4);
 
  // Minimizer
  auto minim = ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"); 

  minim->SetPrintLevel(0);
  minim->SetPrecision(1e-10); 

  // Set minimizer function
  minim->SetFunction(func);

  // Assign variable limits
  minim->SetLimitedVariable(0,"X",parameter[0],0.01,-10,10);
  minim->SetLimitedVariable(1,"Y",parameter[1],0.01,-10,10);
  minim->SetLimitedVariable(2,"Z",parameter[2],0.01,-5,5);
  minim->SetLimitedVariable(3,"T",parameter[3],0.01, 4.7-3.0, 4.7+3.0);

  // Don't draw iterations of minimizer
  flagDraw = 0;

  // Shrink it, babeyyy
  minim->Minimize(); 
  
  // Draw minimal value
  flagDraw = 1;

  // Pull values from minimizer
  const double* x = minim->X();
  cout << "========================================" << endl;
  cout << "\t\tX =" << x[0] << endl;
  cout << "\t\tY =" << x[1] << endl;
  cout << "\t\tZ =" << x[2] << endl;
  cout << "\t\tT =" << x[3] << endl;
  cout << "Minimum: " << devE(x) << endl;
  cout << "========================================" << endl;
}
