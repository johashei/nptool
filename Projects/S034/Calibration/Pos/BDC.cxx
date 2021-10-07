

vector<TVector3*> pos1;
vector<TVector3*> pos2;
vector<TVector3*> posM;

////////////////////////////////////////////////////////////////////////////////
double devR(const double* parameter){
  // Return the standard deviation in Brho
  unsigned int size = pos1.size();

  TVector3 o1(parameter[0],parameter[1],parameter[2]);
  TVector3 o2(parameter[3],parameter[4],parameter[5]);

  double dR;
  TVector3 p1,p2,pB,pM;
  static auto h = new TH1D("h","h", 1000,-0,100);
  h->Reset();
  for(unsigned int i = 0 ; i < size ; i++){
    pM=*(posM[i]);
    p1=*(pos1[i])+o1;
    p2=*(pos2[i])+o2;
    TVector3 BDCDir=p2-p1;
    BDCDir=BDCDir.Unit();
    BDCDir*=(pM.Z()-p2.Z())/BDCDir.Z();
    pB=p2+BDCDir;
    dR=(pB-pM).Mag();
    //cout << (pB-pM).Mag2() << " " << dR << endl;
    double R = sqrt(pB.X()*pB.X()+pB.Y()*pB.Y());
    if(R<15)
      h->Fill(dR); 
  } 
  cout << h->GetMean() << " " << h->GetStdDev() << " " << parameter[0] << endl;
  h->Draw();
  gPad->Update();
  return (h->GetMean()+h->GetStdDev() );
}
////////////////////////////////////////////////////////////////////////////////
void LoadFile(){

  ifstream file("Calibration/Pos/bdc.txt");
  if(!file.is_open()){
    cout << "fail to load file" << endl;
    exit(1);
  }
  else {
    cout <<  "Success opening file" << endl;
  }

  double xm,ym,zm;
  double x1,y1,z1;
  double x2,y2,z2;

  while(file >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> xm >> ym >> zm ){
    auto p1 = new TVector3(x1,y1,z1);
    auto p2 = new TVector3(x2,y2,z2);
    auto pM = new TVector3(xm,ym,zm);
    pos1.push_back(p1);
    pos2.push_back(p2);
    posM.push_back(pM);
  }
  file.close();
}
////////////////////////////////////////////////////////////////////////////////
void BDC(){

  // Data
  LoadFile();

  // Minimizer
  auto min=ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
  auto func=ROOT::Math::Functor(&devR,6); 
  min->SetFunction(func);
  min->SetPrintLevel(0);
  min->SetPrecision(1e-10); 
  double parameter[6]={0,0,-6876.34,0,0,-5876.7};
  devR(parameter);

 // min->SetLimitedVariable(0,"AM",parameter[0],1,-90,90);
  min->SetLimitedVariable(0,"X1",parameter[0],1,-10,10);
  min->SetLimitedVariable(1,"Y1",parameter[1],1,-10,10);
//  min->SetLimitedVariable(2,"Z1",parameter[2],1,parameter[2]-100,parameter[2]+100);
  min->SetFixedVariable(2,"Z1",parameter[2]);
  min->SetLimitedVariable(3,"X2",parameter[3],1,-10,10);
  min->SetLimitedVariable(4,"Y2",parameter[4],1,-10,10);
  //min->SetLimitedVariable(5,"Z2",parameter[5],1,parameter[5]-100,parameter[5]+100);
  min->SetFixedVariable(5,"Z2",parameter[5]);
  min->Minimize(); 

  const double* x = min->X();
  cout << "X1 =" << x[0]<<endl;
  cout << "Y1 =" << x[1]<<endl;
  cout << "Z1 =" << x[2]<<endl;
  cout << "X2 =" << x[3]<<endl;
  cout << "Y2 =" << x[4]<<endl;
  cout << "Z2 =" << x[5]<<endl;
  cout << "Minimum: " << devR(x) << endl;
}
