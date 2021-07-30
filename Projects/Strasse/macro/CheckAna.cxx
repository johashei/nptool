{
gStyle->SetPalette(1);
//TFile *file= new TFile("../../Outputs/Analysis/PhysicsTree.root");
TFile *file= new TFile("../../Outputs/Analysis/configJuly2021_ana.root");
TTree *tree = (TTree*)file->Get("PhysicsTree");

TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
c1->Divide(3,3);
c1->cd(1);
tree->Draw("Ex>>h1(200,-10,10)","","");
c1->cd(2);
tree->Draw("fRC_Vertex_Position_X:fRC_Vertex_Position_Z>>h3(400,-150,250,200,-100,100)","","colz");
c1->cd(3);
tree->Draw("VertexX:VertexZ>>h4(400,-150,250,200,-100,100)","","colz");
c1->cd(4);
tree->Draw("fRC_Vertex_Position_X-VertexX>>h5(100,-2,2)","","");
c1->cd(5);
tree->Draw("fRC_Vertex_Position_Y-VertexY>>h6(100,-2,2)","","");
c1->cd(6);
tree->Draw("fRC_Vertex_Position_Z-VertexZ>>h7(100,-2,2)","","");
c1->cd(7);
//tree->Draw("fRC_Kinetic_Energy[1]-Catana.Energy[0]>>h8(100,-10,10)","Catana.GetEventMultiplicity()==2","");
tree->Draw("fRC_Kinetic_Energy[1]-E1>>h8(100,-50,50)","","");
c1->cd(8);
//tree->Draw("fRC_Kinetic_Energy[1]-Catana.Energy[1]>>h9(100,-10,10)","Catana.GetEventMultiplicity()==2","");
tree->Draw("fRC_Kinetic_Energy[0]-E2>>h9(100,-50,50)","","");
c1->cd(9);
tree->Draw("fRC_Kinetic_Energy[0]:E2>>h10(300,0,300,300,0,300)","","colz");

}
