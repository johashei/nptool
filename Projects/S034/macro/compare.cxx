void compare(){

 auto fz = new TFile("root/zaihong/run0582_RIG20210424_6He.root");
 auto tz = (TTree*) fz->FindObjectAny("rig");
 auto fl = new TFile("root/analysis/test582.root");
 auto tl = (TTree*) fl->FindObjectAny("PhysicsTree");
 tl->AddFriend(tz);

 auto cfdc0= new TCanvas();
 cfdc0->Divide(2,4);
 cfdc0->cd(1);
 tz->Draw("FDC0_X>>h1(1000,-80,80)","FDC0_X>-80 && FDC0_X<80");
 cfdc0->cd(2);
 tz->Draw("FDC0_Y>>h2(1000,-80,80)","FDC0_X>-80 && FDC0_Y<80");
 cfdc0->cd(3);
 tl->Draw("SamuraiFDC0.PosX>>h3(1000,-80,80)","SamuraiFDC0.PosX>-80 && SamuraiFDC0.PosX<80");
 cfdc0->cd(4);
 tl->Draw("SamuraiFDC0.PosY>>h4(1000,-80,80)","SamuraiFDC0.PosY>-80 && SamuraiFDC0.PosY<80");
 cfdc0->cd(5);
 tl->Draw("SamuraiFDC0.PosX:FDC0_X>>h5(1000,-80,80,1000,-80,80)","SamuraiFDC0.PosX>-80 && SamuraiFDC0.PosX<80","colz");
 auto l = new TLine(-80,-80,80,80);
 l->SetLineColor(kOrange+7);
 l->SetLineWidth(2);
 l->Draw();
 cfdc0->cd(6);
 tl->Draw("SamuraiFDC0.PosY:FDC0_Y>>h6(1000,-80,80,1000,-80,80)","SamuraiFDC0.PosY>-80 && SamuraiFDC0.PosY<80","colz");
 l->Draw();

 cfdc0->cd(7);
 tl->Draw("SamuraiFDC0.PosX-FDC0_X>>h7(1000,-5,5)","SamuraiFDC0.PosX>-80 && SamuraiFDC0.PosX<80","");
 cfdc0->cd(8);
 tl->Draw("SamuraiFDC0.PosY-FDC0_Y>>h8(1000,-5,5)","SamuraiFDC0.PosY>-80 && SamuraiFDC0.PosY<80","");
 
 auto cfdc2= new TCanvas();
 cfdc2->Divide(2,4);
 cfdc2->cd(1);
 tz->Draw("FDC2_X>>g1(1000,-2000,2000)","FDC2_X>-2000 && FDC2_X<2000");
 cfdc2->cd(2);
 tz->Draw("FDC2_Y>>g2(1000,-2000,2000)","FDC2_X>-2000 && FDC2_Y<2000");
 cfdc2->cd(3);
 tl->Draw("SamuraiFDC2.PosX>>g3(1000,-2000,2000)","SamuraiFDC2.PosX>-2000 && SamuraiFDC2.PosX<2000&&SamuraiFDC2.devX<10");
 cfdc2->cd(4);
 tl->Draw("SamuraiFDC2.PosY>>g4(1000,-2000,2000)","SamuraiFDC2.PosY>-2000 && SamuraiFDC2.PosY<2000&&SamuraiFDC2.devY<10");
 cfdc2->cd(5);
 tl->Draw("SamuraiFDC2.PosX:FDC2_X>>g5(1000,-2000,2000,1000,-2000,2000)","SamuraiFDC2.PosX>-2000 && SamuraiFDC2.PosX<2000&&SamuraiFDC2.devX<10","colz");
 l = new TLine(-2000,-2000,2000,2000);
 l->SetLineColor(kOrange+7);
 l->SetLineWidth(2);
 l->Draw();

 cfdc2->cd(6);
 tl->Draw("SamuraiFDC2.PosY:FDC2_Y>>g6(1000,-2000,2000,1000,-2000,2000)","SamuraiFDC2.PosY>-2000 && SamuraiFDC2.PosY<2000&&SamuraiFDC2.devX<10","colz");
 l->Draw();

 cfdc2->cd(7);
 tl->Draw("SamuraiFDC2.PosX-FDC2_X>>g7(1000,-5,5)","SamuraiFDC2.PosX>-2000 && SamuraiFDC2.PosX<2000&&SamuraiFDC2.devX<10","");
 cfdc2->cd(8);
 tl->Draw("SamuraiFDC2.PosY-FDC2_Y>>g8(1000,-5,5)","SamuraiFDC2.PosY>-2000 && SamuraiFDC2.PosY<2000&&SamuraiFDC2.devY<10","");
 
 auto cfdc2b= new TCanvas();
 cfdc2b->Divide(2,4);
 cfdc2b->cd(1);
 tz->Draw("FDC2_ThetaX>>i1(1000,-2,2)","FDC2_ThetaX>-1000");
 cfdc2b->cd(2);
 tz->Draw("FDC2_ThetaY>>i2(1000,-2,2)","FDC2_ThetaY>-1000");
 cfdc2b->cd(3);
 tl->Draw("SamuraiFDC2.ThetaX>>i3(1000,-2,2)","SamuraiFDC2.ThetaX>-10000");
 cfdc2b->cd(4);
 tl->Draw("SamuraiFDC2.PhiY>>i4(1000,-2,2)","SamuraiFDC2.PhiY>-1000");
 cfdc2b->cd(5);
 tl->Draw("SamuraiFDC2.ThetaX:FDC2_ThetaX>>i5(1000,-2,2,1000,-2,2)","FDC2_ThetaX>-1000&&SamuraiFDC2.ThetaX>-10000","colz");
 l = new TLine(-2,-2,2,2);
 l->SetLineColor(kOrange+7);
 l->SetLineWidth(2);
 l->Draw();

 cfdc2b->cd(6);
 tl->Draw("SamuraiFDC2.PhiY:FDC2_ThetaY>>i6(1000,-2,2,1000,-2,2)","FDC2_ThetaY>-1000&&SamuraiFDC2.PhiY>-1000","colz");
 l->Draw();

 cfdc2b->cd(7);
 tl->Draw("SamuraiFDC2.ThetaX-FDC2_ThetaX>>i7(1000,-0.1,0.1)","FDC2_ThetaX>-1000&&SamuraiFDC2.ThetaX>-10000","");
 cfdc2b->cd(8);
 tl->Draw("SamuraiFDC2.PhiY-FDC2_ThetaY>>i8(1000,-0.1,0.1)","FDC2_ThetaY>-1000&&SamuraiFDC2.PhiY>-1000","");

 auto cbroh= new TCanvas();
 cbroh->Divide(2,2);
 cbroh->cd(1);
 tl->Draw("Brho>>b1(1000,0,8)","Brho>0 && SamuraiFDC2.devX<10 && SamuraiFDC2.devX<10");
 cbroh->cd(2);
 tz->Draw("rig>>b2(1000,0,8000)","rig>0");
 cbroh->cd(3);
 tl->Draw("Brho:rig/1000.>>b3(1000,0,8,1000,0,8)","rig>0 && Brho>0","colz");
 l = new TLine(0,0,8,8);
 l->SetLineColor(kOrange+7);
 l->SetLineWidth(2);
 l->Draw();
 cbroh->cd(4);
 tz->Draw("rig/1000.:FragID>>b4(1000,-2000,2000,31,-1,30)","","colz");

}
