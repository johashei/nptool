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

    //Define normal vector for the MG# of detection
    //[[[[ UPDATE WITH NEW MG POSITIONS FROM SURVEY ]]]]
    switch(detnum[i]){
      case 1:
        MugastNormal.SetXYZ(-0.453915, +0.455463, -0.765842);
        break;
      case 2:
	MugastNormal.SetXYZ(-0.642828, +0.000000, -0.766010);
        break;
      case 3:
	MugastNormal.SetXYZ(-0.454594, -0.450670, -0.768271);
        break;
      case 4:
	MugastNormal.SetXYZ(-0.002437, -0.638751, -0.769409);
        break;
      case 5:
	MugastNormal.SetXYZ(+0.452429, -0.454575, -0.767248);
        break;
      case 7:
	MugastNormal.SetXYZ(+0.443072, +0.443265, -0.779232);
        break;
      default:
	cout << "ERROR:: Invalid DetNum " << detnum[i] << " at event " << i << endl;
        return 1; // Exit code
    }

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
        cout << "ERROR:: Invalid DetNum " << detnum[i] << " at event " << i << endl;
        return 1; // Exit code
    }
  }
  //End loop over events

  //Select MG# being minimized
  double mean = 0;
  double stddev = 0;
  switch(mgSelect){
    case 0:
      mean   = h->GetMean(); 
      stddev = h->GetStdDev(); 
      break;
    case 1:
      mean   = h1->GetMean(); 
      stddev = h1->GetStdDev(); 
      break;
    case 2:
      mean   = h2->GetMean(); 
      stddev = h2->GetStdDev(); 
      break;
    case 3:
      mean   = h3->GetMean(); 
      stddev = h3->GetStdDev(); 
      break;
    case 4:
      mean   = h4->GetMean(); 
      stddev = h4->GetStdDev(); 
      break;
    case 5:
      mean   = h5->GetMean(); 
      stddev = h5->GetStdDev(); 
      break;
    case 6:
      mean   = h7->GetMean(); 
      stddev = h7->GetStdDev(); 
      break;
    default:
      cout << "ERROR:: Invalid MG# selection! -> " << mgSelect << endl;
      return 1; // Exit code
  }

  //Write vals to screen
  cout << "Mean: " << mean 
       << "    StdDev: " << stddev 
       << "    Thickness: " << parameter[3] << " um" 
       << endl;

  //Draw histogram(s)
  h->Draw();
  if(flagDraw){ InitiliseCanvas(); }

  /*
  cout << pow(h->GetMean()-refE,2)  << " + " 
       << pow(0.1*h->GetStdDev(),2) << " = " 
       << pow(h->GetMean()-refE,2) + pow(0.1*h->GetStdDev(),2) 
       << endl;
  */

  //Adapt the metric as needed
  return sqrt( pow(mean-refE,2) + pow(0.1*stddev,2) );
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
 
  //Minimizer
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

  //Shrink it, babeyyy
  minim->Minimize(); 
  
  //Draw minimal value
  flagDraw = 1;

  //Pull values from minimizer
  const double* x = minim->X();
  cout << "==================================================" << endl;
  cout << "=---------------- FINAL PEAK FITS ---------------=" << endl;
  cout << "==================================================" << endl;
    devE(x);
  cout << "==================================================" << endl;
  cout << "=------------ RESULTS OF MINIMIZATION -----------=" << endl;
  cout << "==================================================" << endl;
  cout << "\t\tX =" << x[0] << endl;
  cout << "\t\tY =" << x[1] << endl;
  cout << "\t\tZ =" << x[2] << endl;
  cout << "\t\tT =" << x[3] << endl;
  cout << "==================================================" << endl;
}
