void CheckMinos(){

  auto f1 = new TFile("root/analysis/run582CheckMinos.root","READ");
  auto tree1 = (TTree*) f1->FindObjectAny("PhysicsTree");
  auto f2 = new TFile("root/analysis/run824CheckMinos.root","READ");
  auto tree2 = (TTree*) f2->FindObjectAny("PhysicsTree");

  auto c = new TCanvas();
  c->Divide(2,2);
  c->cd(1);
  tree1->Draw("Minos.Z_Vertex>>h1(500,-150,250)","sqrt(X_Vertex*X_Vertex+Y_Vertex*Y_Vertex)<15");
  c->cd(2);
  tree1->Draw("Minos.Delta_Vertex>>h2(100,0,10)");
  c->cd(3);
  tree2->Draw("Minos.Z_Vertex>>h3(500,-150,250)","sqrt(X_Vertex*X_Vertex+Y_Vertex*Y_Vertex)<15");
  c->cd(4);
  tree2->Draw("Minos.Delta_Vertex>>h4(100,0,10)");




}
