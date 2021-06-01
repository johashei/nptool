

vector<TVector3*> pos0;
vector<TVector3*> dir0;
vector<TVector3*> pos2;
vector<TVector3*> dir2;
SamuraiFieldMap field;

////////////////////////////////////////////////////////////////////////////////
double devB(const double* parameter){
  // Return the standard deviation in Brho
  unsigned int size = pos0.size();
  TVector3 o0(parameter[0],parameter[1],parameter[2]);
  TVector3 o2(parameter[3],parameter[4],parameter[5]);
  TVector3 p0,p2;
  field.SetFDC2R(parameter[5]);
  double Brho;
  static auto h = new TH1D("h","h", 1000,0,10);
  h->Reset();
  for(unsigned int i = 0 ; i < 1000 ; i++){
    p0=*(pos0[i])+o0;
    p2=*(pos2[i])+o2;

    if(dir0[i]->Z()>0.8)
    h->Fill(field.FindBrho(p0,*dir0[i],p2,*dir2[i])); 
  } 
  o0.Print();
  o2.Print();
  cout << h->GetStdDev() << endl;
  return h->GetStdDev();
}
////////////////////////////////////////////////////////////////////////////////
void LoadFile(){

  ifstream file("Calibration/Pos/fdc.txt");
  if(!file.is_open()){
    cout << "fail to load file" << endl;
    exit(1);
  }
  else {
    cout <<  "Success opening file" << endl;
  }

  double x0,y0,z0;
  double dx0,dy0,dz0;
  double x2,y2,z2;
  double dx2,dy2,dz2;

  while(file >> x0 >> y0 >> z0 >> dx0 >> dy0 >> dz0 >> x2 >> y2 >> z2 >> dx2 >> dy2 >> dz2){
    auto p0 = new TVector3(x0,y0,z0);
    auto d0 = new TVector3(dx0,dy0,dz0);
    auto p2 = new TVector3(x2,y2,z2);
    auto d2 = new TVector3(dx2,dy2,dz2);
    pos0.push_back(p0);
    dir0.push_back(d0);
    pos2.push_back(p2);
    dir2.push_back(d2);
  }
  cout << " Val " << pos0.size() << " Value Loaded" << endl;
  file.close();
}
////////////////////////////////////////////////////////////////////////////////
void FDC(){

  // Data
  LoadFile();

  // Field map
  field.LoadMap(30*deg,"field_map/180702-2,40T-3000.table.bin",10);
  field.SetFDC2Angle((59.930-90.0)*deg);

  double parameter[6]={0,0,-3456.52,-252.55,0,4123.47};
  cout << devB(parameter) << endl;

  // Minimizer
  auto min=ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
  auto func=ROOT::Math::Functor(&devB,6); 
  min->SetFunction(func);
  min->SetPrintLevel(0);
  min->SetPrecision(1e-3); 
  min->SetLimitedVariable(0,"X0",0,1,-10,10);
  min->SetLimitedVariable(1,"Y0",0,1,-10,10);
  min->SetLimitedVariable(2,"Z0",-3372.07,1,-3500,-3400);
  min->SetLimitedVariable(3,"X2",-252.55,1,-260,-240);
  min->SetLimitedVariable(4,"Y2",0,1,-10,10);
  min->SetLimitedVariable(5,"Z2",4123.47,1,4000,4200);
  min->Minimize(); 
  const double* x = min->X();
  cout << "X0 =" << x[0]<<endl;
  cout << "Y0 =" << x[1]<<endl;
  cout << "Z0 =" << x[2]<<endl;
  cout << "X2 =" << x[3]<<endl;
  cout << "Y2 =" << x[4]<<endl;
  cout << "Z2 =" << x[5]<<endl;

}
