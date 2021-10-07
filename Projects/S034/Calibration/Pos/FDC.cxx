

vector<TVector3*> pos0;
vector<TVector3*> posM;
vector<TVector3*> pos2;
vector<TVector3*> dir2;
SamuraiFieldMap field;

////////////////////////////////////////////////////////////////////////////////
double devB(const double* parameter){
  // Return the standard deviation in Brho
  unsigned int size = pos0.size();
  TVector3 oM(parameter[1],parameter[2],parameter[3]);
  TVector3 o0(parameter[4],parameter[5],parameter[6]);
  TVector3 o2(parameter[7],parameter[8],parameter[9]);
  TVector3 pM,p0,d0,p2;
  oM.Print();
  o0.Print();
  o2.Print();
  field.SetFDC2R(parameter[6]);
  static auto h = new TH1D("h","h", 1000,0,10);
  h->Reset();
  for(unsigned int i = 0 ; i < size ; i++){
    pM=*(posM[i]);
    pM.RotateZ(parameter[0]*M_PI/180.);
    pM+=oM;
    p0=*(pos0[i])+o0;
    p2=*(pos2[i])+o2;
    d0=(p0-pM).Unit();
//    p2.Print();
    if(d0.Z()>0.9)
      h->Fill(field.FindBrho(p0,d0,p2,*dir2[i])); 
  } 
  cout << h->GetStdDev() << " " << parameter[0] << endl;
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

  double xm,ym,zm;
  double x0,y0,z0;
  double x2,y2,z2;
  double dx2,dy2,dz2;

  while(file >> x0 >> y0 >> z0 >> xm >> ym >> zm >> x2 >> y2 >> z2 >> dx2 >> dy2 >> dz2){
    auto p0 = new TVector3(x0,y0,z0);
    auto pM = new TVector3(xm,ym,zm);
    auto p2 = new TVector3(x2,y2,z2);
    auto d2 = new TVector3(dx2,dy2,dz2);
    pos0.push_back(p0);
    posM.push_back(pM);
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

  // Minimizer
  auto min=ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
  auto func=ROOT::Math::Functor(&devB,10); 
  min->SetFunction(func);
  min->SetPrintLevel(0);
  min->SetPrecision(1e-6); 
  double parameter[10]={40.6,0,0,-4657.39, 0, 0,-3372.07,-252.55,0,4123.47};
  devB(parameter);

 // min->SetLimitedVariable(0,"AM",parameter[0],1,-90,90);
  min->SetFixedVariable(0,"AM",parameter[0]);
  min->SetLimitedVariable(1,"XM",parameter[1],1,-10,10);
  min->SetLimitedVariable(2,"YM",parameter[2],1,-10,10);
  min->SetLimitedVariable(3,"ZM",parameter[3],1,-4650,-4660);
  min->SetLimitedVariable(4,"X0",parameter[4],1,-10,10);
  min->SetLimitedVariable(5,"Y0",parameter[5],1,-10,10);
  min->SetLimitedVariable(6,"Z0",parameter[6],1,-3370,-3375);
  min->SetLimitedVariable(7,"X2",parameter[7],1,-260,-240);
  min->SetLimitedVariable(8,"Y2",parameter[8],1,-10,10);
  min->SetLimitedVariable(9,"Z2",parameter[9],1,4120,4125);
  /*min->SetFixedVariable(1,"XM",parameter[1]);
  min->SetFixedVariable(2,"YM",parameter[2]);
  min->SetFixedVariable(3,"ZM",parameter[3]);
  min->SetFixedVariable(4,"X0",parameter[4]);
  min->SetFixedVariable(5,"Y0",parameter[5]);
  min->SetFixedVariable(6,"Z0",parameter[6]);
  min->SetFixedVariable(7,"X2",parameter[7]);
  min->SetFixedVariable(8,"Y2",parameter[8]);
  min->SetFixedVariable(9,"Z2",parameter[9]);
 */
  min->Minimize(); 
  const double* x = min->X();
  cout << "AM =" << x[0]<<endl;
  cout << "XM =" << x[1]<<endl;
  cout << "YM =" << x[2]<<endl;
  cout << "ZM =" << x[3]<<endl;
  cout << "X0 =" << x[4]<<endl;
  cout << "Y0 =" << x[5]<<endl;
  cout << "Z0 =" << x[6]<<endl;
  cout << "X2 =" << x[7]<<endl;
  cout << "Y2 =" << x[8]<<endl;
  cout << "Z2 =" << x[9]<<endl;
  cout << "Minimum: " << devB(x) << endl;
}
