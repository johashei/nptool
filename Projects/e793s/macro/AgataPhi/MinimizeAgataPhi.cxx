#include "MinimizeAgataPhi.h"
////////////////////////////////////////////////////////////////////////////////
double devE(const double * parameter) {

  //Other variable initilizations
  unsigned int size = GammaE.size();
  double rotationx = parameter[0];
  double rotationy = parameter[1];
  double rotationz = parameter[2];

  //Initilize histogram
  h->Reset();

  //Now that initial range is wide, crop to single peak
  //if(refE==0.143){
  //  h->SetAxisRange(-1.0, +1.0, "X");
  //}

  //Initilize results array
  //    7 => Sum in 0 and them MG's in 1-6
  //    5 => Mean, MeanErr, StdDev, StdDevErr, Chi2/NDF
  double FitResultMatrix[5];

  //Loop over events
  for (unsigned int i = 0; i < size; i++) {
    //Construct LorentzV
    TLorentzVector GammaLV;
    GammaLV.SetPx(GammaP[i]->X());
    GammaLV.SetPy(GammaP[i]->Y());
    GammaLV.SetPz(GammaP[i]->Z());
    GammaLV.SetE(GammaE[i]);
    //Construct beta vetoc, with rotation
    TVector3 beta(Beta[i]->X(), Beta[i]->Y(), Beta[i]->Z());
    beta.RotateX(rotationx);
    beta.RotateY(rotationy);
    beta.RotateZ(rotationz);
    //Boost & fill
    GammaLV.Boost(beta);
    h->Fill(GammaLV.Energy());
  }

  //Initilise, Draw & Fit histograms
  InitiliseCanvas(FitResultMatrix);

  //Write vals to screen
  if (flagDraw) {cout << "==================================================" << endl;}

  /* Metric */
  double multiplier = 100.; //0.08;
  double metric = FitResultMatrix[2];

  WriteToCout(FitResultMatrix, parameter[3], metric);
  if (flagDraw) { WriteToFile(FitResultMatrix, parameter, metric); }
  
  return metric;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void MinimizeAgataPhi() {

  ////filename = "XYZE_ShiftMG3_GateOn0143.txt";
  filename = "RotBeta_17Feb_Large.txt";
  //refE = 0.1427; // the energy of the selected states
  //refE = 0.279; // the energy of the selected states
  //refE = 0.449; // the energy of the selected states
  refE = 1.838; // the energy of the selected states

  //Read data in
  LoadFile();

  //Output formatting
  cout << fixed << showpoint << setprecision(6) << showpos;

  //Start timer
  auto start = high_resolution_clock::now();

  //Start with pi rotation in y
  double parameter[3] = {0.0,M_PI,0.0};

  //Don't draw iterations of minimizer
  flagDraw = 0;
  gROOT -> SetBatch(kTRUE);

  //Initial pass through
  devE(parameter);

  //Function with 3 parameters, ROtX, RotY, RotZ
  auto func = ROOT::Math::Functor( & devE, 3);

  //Initilise minimizer
  auto minim = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  minim -> SetMaxFunctionCalls(100000000); // used by Minuit and Minuit2 
  minim -> SetMaxIterations(100000000); // used by GSL
  minim -> SetPrintLevel(3);
  minim -> SetPrecision(1e-18);

  //Set minimizer function
  minim -> SetFunction(func);

  //Assign variable limits
  minim -> SetLimitedVariable(0, "RotationX", parameter[0], 0.001,-0.1*M_PI,+0.1*M_PI);
  //minim -> SetFixedVariable(0, "RotationX", parameter[0]);
  minim -> SetLimitedVariable(0, "RotationY", parameter[1], 0.001,+0.9*M_PI,+1.1*M_PI);
  //minim -> SetFixedVariable(1, "RotationY", parameter[1]);
  minim -> SetLimitedVariable(2, "RotationZ", parameter[2], 0.001,-0.1*M_PI,+0.1*M_PI);
  //minim -> SetFixedVariable(2, "RotationZ", parameter[2]);

  //Don't draw iterations of minimizer
  flagDraw = 0;
  gROOT -> SetBatch(kTRUE);

  // Rebin for high energy gammas
  if(refE>1.0){ h->Rebin(2);}

  //Shrink it, babeyyy
  minim -> Minimize();

  //Draw minimal value
  flagDraw = 1;
  gROOT -> SetBatch(kFALSE);

  //Pull values from minimizer
  const double * x = minim -> X();
  cout << "==================================================" << endl;
  cout << "=---------------- FINAL PEAK FITS ---------------=" << endl;
  cout << "==================================================" << endl;
    devE(x);
  cout << "==================================================" << endl;
  cout << "=------------ RESULTS OF MINIMIZATION -----------=" << endl;
  cout << "\t\tRotateX = " << x[0] << " " << endl;
  cout << "\t\tRotateY = " << x[1] << " " << endl;
  cout << "\t\tRotateZ = " << x[2] << " " << endl;
  cout << "==================================================" << endl;

  //Stop timer
  auto stop = high_resolution_clock::now();

  //Close output file
  file.close();

  //Program runtime
  auto duration = duration_cast < seconds > (stop - start);
  cout << " Runtime = " << duration.count() << " s" << endl;
}
