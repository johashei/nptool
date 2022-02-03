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
  auto Ex1p=new TH1D("ex1p","ex1p",300,-10,40); 
  PS1p.SetDecay(W,3,masses1p);

  for(unsigned int i = 0 ; i < 1000000 ; i++){
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

  Ex1p->Draw("");
}

