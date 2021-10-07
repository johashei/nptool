void IDz(){
  auto fz = new TFile("root/zaihong/run0582_RIG20210424_6He.root");
  auto tz = (TTree*) fz->FindObjectAny("rig");
 
  tz->Draw("rig/1000.>>h(1000,0,8)","","");  
  tz->Draw("rig/1000.>>hn(1000,0,8)","FragID==-1","same");  
  tz->Draw("rig/1000.>>hd(1000,0,8)","FragID==12","same");  
  tz->Draw("rig/1000.>>ht(1000,0,8)","FragID==13","same");  
  tz->Draw("rig/1000.>>h4He(1000,0,8)","FragID==24","same");  
  tz->Draw("rig/1000.>>h6He(1000,0,8)","FragID==26","same");  
  auto hn = (TH1*) gDirectory->FindObjectAny("hn");
  auto hd = (TH1*) gDirectory->FindObjectAny("hd");
  auto ht = (TH1*) gDirectory->FindObjectAny("ht");
  auto h4He = (TH1*) gDirectory->FindObjectAny("h4He");
  auto h6He = (TH1*) gDirectory->FindObjectAny("h6He");

  hn->SetFillStyle(1001);
  hd ->SetFillStyle(1001);
  ht  ->SetFillStyle(1001);
  h4He ->SetFillStyle(1001);
  h6He ->SetFillStyle(1001);

  hn->SetFillColor(kBlack);
  hd ->SetFillColor(kAzure+7);
  ht  ->SetFillColor(kOrange+7);
  h4He ->SetFillColor(kGreen-3);
  h6He ->SetFillColor(kMagenta-3);
  hn->SetLineColor(kBlack);
  hd ->SetLineColor(kAzure+7);
  ht  ->SetLineColor(kOrange+7);
  h4He ->SetLineColor(kGreen-3);
  h6He ->SetLineColor(kMagenta-3);
  gPad->Update();

}
