Double_t f_semi(Double_t *x, Double_t *par){
  Float_t xx = x[0];
  Double_t f;
  Double_t a = TMath::Power(par[0],2);
  Double_t b = TMath::Power(xx-par[1],2);
  if(a > b){ f = par[2]*TMath::Sqrt(a-b); }
  else{ f = 0; }
  return f;
}


void ThreeBodyBreakup(){
  TGenPhaseSpace PS1p;

  //47K(d,p)47K+p+n
  NPL::Reaction R("47K(d,p)48K@355");
  NPL::Particle K47("47K");
  NPL::Particle p("1H");
  NPL::Particle d("2H");
  NPL::Particle n("n");

  TLorentzVector target(0,0,0,d.Mass());
  K47.SetKineticEnergy(354.75);
  double E = K47.GetGamma()*K47.Mass();
  double pc=sqrt(E*E-K47.Mass()*K47.Mass());
  TLorentzVector beam(0,0,pc,K47.GetGamma()*K47.Mass());
  TLorentzVector W = target+beam;

  double masses1p[3]={p.Mass(),K47.Mass(),n.Mass()};
  //cout << p.Mass() << " " << K47.Mass()<< " " << n.Mass() << endl;
  auto Ex1p=new TH1D("ex1p","ex1p",500,-10,40); 
  PS1p.SetDecay(W,3,masses1p);

  for(unsigned int i = 0 ; i < 10000000 ; i++){
     // 1p case
     double weight = PS1p.Generate();
     auto pd = PS1p.GetDecay(0);
     double theta = pd->Vect().Angle(TVector3(0,0,1));
     double beta  = pd->Beta();
     p.SetBeta(beta);
     double ELab=p.GetEnergy();
     double ex=R.ReconstructRelativistic(ELab,theta);
     //cout << ELab << " " << ex << endl;
     Ex1p->Fill(ex,weight); 
    }

  TCanvas* c_Ex1p = new TCanvas("c_Ex1p","c_Ex1p",1000,500);
  Ex1p->Draw();

//  TCanvas* c_semi = new TCanvas("c_semi","c_semi",1000,500);
  TF1 *f1 = new TF1("f_semi",f_semi,-10,40,3);
  f1->SetParameters(6,10,100);
  f1->SetParNames("Radius","Centre","VertScale");

  TFitResultPtr result = Ex1p->Fit("f_semi","S+");

  f1->Draw("same");

  cout << "Chi2 = " << result->Chi2() 
       << "  NDF  = " << result->Ndf() 
       << "  Chi2/NDF = " << result->Chi2()/result->Ndf() 
       << endl;

  cout << "==========================================" << endl;
  cout << " Three-Body Breakup, B(x);" << endl;
  cout << "     B(x) = " 
       << result->Value(2) << " * SQRT( (" 
       << result->Value(0) << ")^2 - (x - " 
       << result->Value(1) << ")^2 )" 
       << endl;

}

