
double ref = 0.143; // the energy of the selected states
vector<TVector3*> pos;
vector<double> energy;
NPL::EnergyLoss CD2("proton_CD2.G4table","G4Table",100);
NPL::EnergyLoss Al("proton_Al.G4table","G4Table",100);

////////////////////////////////////////////////////////////////////////////////
double devE(const double* parameter){
  
  // My reaction at 8A MeV*47A=376 MeV
  static NPL::Reaction reaction("47K(d,p)48K@376");
  // Return the standard deviation in Brho
  unsigned int size = pos1.size();

  TVector3 offset(parameter[0],parameter[1],parameter[2]);

  double dE,Theta;
  TVector3 dir;
  static auto h = new TH1D("h","h", 1000,-0,100);
  h->Reset();
  for(unsigned int i = 0 ; i < size ; i++){
    dir=*(pos[i])-o1;
    double Theta= dir.Angle(TVector3(0,0,1));
    double Energy = energy[i];
    // do some energy loss correction here:
    // This need to take parameter[4]*0.5*micrometer as the target thickness
    Energy=CD2.EvaluateInitialEnergy(Energy,0,5*parameter[4]*micrometer,Theta);

    double Ex=reaction->ReconstructRelativistic(Energy,Theta);
    //cout << (pB-pM).Mag2() << " " << dR << endl;
      h->Fill(Ex-ref); 
  } 
  cout << h->GetMean() << " " << h->GetStdDev() << " " << parameter[4] << endl;
  h->Draw();
  gPad->Update();
  // adapt the metric as needed
  return (h->GetMean()+h->GetStdDev() );
}
////////////////////////////////////////////////////////////////////////////////
void LoadFile(){

  ifstream file("myFile.txt");
  if(!file.is_open()){
    cout << "fail to load file" << endl;
    exit(1);
  }
  else {
    cout <<  "Success opening file" << endl;
  }

  double x,y,z,e;

  while(file >> x >> y >> z >> e ){
    auto p1 = new TVector3(x,y,z);
    double e;
    pos.push_back(p1);
    energy.push_back(e);
  }
  file.close();
}
////////////////////////////////////////////////////////////////////////////////
void BDC(){

  // Data
  LoadFile();

  // Minimizer
  auto min=ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
  // Function with 4 parameter XYZ and Target thickness
  auto func=ROOT::Math::Functor(&devE,4); 
  min->SetFunction(func);
  min->SetPrintLevel(0);
  min->SetPrecision(1e-10); 
  // Start value: Center beam (0,0,0) and 4.7um target 0.5mg/cm2 target
  double parameter[6]={0,0,0,4.7};
  devR(parameter);

 // min->SetLimitedVariable(0,"AM",parameter[0],1,-90,90);
  min->SetLimitedVariable(0,"X",parameter[0],0.1,-10,10);
  min->SetLimitedVariable(1,"Y",parameter[1],0.1,-10,10);
  min->SetLimitedVariable(2,"Z",parameter[2],0.1,-5,5);
  min->SetLimitedVariable(3,"T",parameter[3],0.1,parameter[2]-0.1*parameter[2],parameter[2]+0.1*parameter[2]);
  min->Minimize(); 

  const double* x = min->X();
  cout << "X =" << x[0]<<endl;
  cout << "Y =" << x[1]<<endl;
  cout << "Z =" << x[2]<<endl;
  cout << "Y =" << x[3]<<endl;
  cout << "Minimum: " << devR(x) << endl;
}
