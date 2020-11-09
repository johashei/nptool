#include"NPDCReconstruction.h"
#include"NPSystemOfUnits.h"
using namespace NPUNITS;
R__LOAD_LIBRARY(libNPTrackReconstruction.dylib)
  //
  void DrawTrack(double a, double b, int style);
  //
  void DrawCircle(double x, double z, double r);
  //
  void DCTest(){
    // load the lib
    auto c = new TCanvas();
    c->Divide(2,1);
    TRandom3 rand;
    // The ronconstruction utility object
    NPL::DCReconstruction dcr;
    vector<double> x,z,r;
    // maximum drift length
    double maxR = 10; 
    double stepZ = 50;

    // a drawing artefact
    // give axis and background to ResolvePlane
    auto bck = new TH2F("bck","bck",10,-150,150,10,-150,150);

    // infinite loop
    while(true){
      // Check 2D minimisation
      c->cd(1);
      x.clear();z.clear();r.clear();
      // create a random trajectory in 2D //
      // a and b are the parameter of the 2D track
      double a = rand.Uniform(-10,10);
      double b = rand.Uniform(-100,100); 
      DrawTrack(a,b,2); 
      // Generate n points from the line
      int n = (int) rand.Uniform(3,20);
      for(unsigned int i = 0 ; i < n ; i++){
        double Z = i*stepZ-100;
        z.push_back(Z); 
        double R = rand.Uniform(-maxR,maxR);
        r.push_back(R);
        double X = (Z-b)/a+R;
        x.push_back(X);
        DrawCircle(X,Z,R);
      }
      // add second track for pile up
      a = rand.Uniform(-10,10);
      b = rand.Uniform(-100,100); 
      DrawTrack(a,b,3); 
      // Generate n points from the line
      n = (int) rand.Uniform(3,20);
      for(unsigned int i = 0 ; i < n ; i++){
        double Z = i*stepZ-100;
        z.push_back(Z); 
        double R = rand.Uniform(-maxR,maxR);
        r.push_back(R);
        double X = (Z-b)/a+R;
        x.push_back(X);
        DrawCircle(X,Z,R);
      }

      double x0,x100,af,bf;
      double res = dcr.BuildTrack2D(x,z,r,x0,x100,af,bf);
      cout << res << endl;
      DrawTrack(af,bf,1);

      // Check 2D track crossing point
      c->cd(2);
      bck->Draw();
      double ThetaL=-60*deg;
      double ThetaH= 60*deg;
      TVector3 L(rand.Uniform(-100,100),0,0);
      TVector3 H(rand.Uniform(-100,100),0,0);
      L.RotateZ(ThetaL);
      H.RotateZ(ThetaH);
      // direction of U and V wire
      TVector3 u(0,1,0);
      u.RotateZ(ThetaL);

      TVector3 v(0,1,0);
      v.RotateZ(ThetaH);

      // Compute the coeff of the two line of vecotr u (v) going through H (L)
      // dv : y = av*x+bv
      double av = v.Y()/v.X();
      double bv = H.Y() - av*H.X();

      DrawTrack(av,bv,1);
      // du : y = au*x+bu
      double au = u.Y()/u.X();
      double bu = L.Y() - au*L.X();

      DrawTrack(au,bu,1);
      TVector3 P;
      dcr.ResolvePlane(L,ThetaL,H,ThetaH,P);
      auto Lmk = new TMarker(L.X(),L.Y(),23);
      auto Hmk = new TMarker(H.X(),H.Y(),27);
      auto Pmk = new TMarker(P.X(),P.Y(),29);
      Lmk->Draw();
      Hmk->Draw();
      Pmk->Draw();

      c->Update();
      gPad->WaitPrimitive();
    }
  }
// draw a track based on z=ax+b
void DrawTrack(double a, double b, int style){
  double x[2]={-100,100};
  double z[2]={a*-100+b,a*100+b};
  auto l = new TGraph(2,x,z);
  l->SetLineColor(kBlack);
  if(style>1)
    l->SetLineWidth(2);
  l->SetLineStyle(style);
  if(style==2)
    l->Draw("al");
  else
    l->Draw("l");
}
// draw a circle
void DrawCircle(double x, double z, double r){
  auto e = new TEllipse(x,z,r,r);
  e->SetLineColor(kBlack);
  e->Draw("same"); 
}
