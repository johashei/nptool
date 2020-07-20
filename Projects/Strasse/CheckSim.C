{
gStyle->SetPalette(1);
//TFile *file= new TFile("../../Outputs/Simulation/test_ref.root");
TFile *file= new TFile("../../Outputs/Simulation/test_optimized.root");
TTree *tree = (TTree*)file->Get("SimulatedTree");

TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
c1->Divide(4,4);
c1->cd(1);
tree->Draw("fRC_Beam_Reaction_Energy:fRC_Vertex_Position_Z>>h1(200,-3.1,3.1,220,,)","","colz");
c1->cd(2);
tree->Draw("fRC_Vertex_Position_Y:fRC_Vertex_Position_Z>>h2","","");
c1->cd(3);
tree->Draw("fRC_Vertex_Position_X:fRC_Vertex_Position_Z>>h3","","");
c1->cd(4);
tree->Draw("fDetected_Position_Y:fDetected_Position_X>>h4","","colz");
c1->cd(5);
tree->Draw("fDetected_Position_Y:fDetected_Position_Z>>h5(1200,-70,230,280,-70,70)","","colz");
c1->cd(6);
tree->Draw("fInner_TE_StripNbr:fDetected_Position_Z>>h6","","colz");
c1->cd(7);
tree->Draw("fInner_LE_StripNbr:fDetected_Position_X>>h7","","colz");
c1->cd(8);
tree->Draw("Strasse.GetOuterMultTEnergy()+Strasse.GetInnerMultTEnergy()>>h8(6,0,6)","","");
c1->cd(9);
tree->Draw("fInner_LE_StripNbr>>h9(1250,0,1250)","","");
c1->cd(10);
tree->Draw("fInner_TE_StripNbr>>h10(1250,0,1250)","","");
c1->cd(11);
tree->Draw("fOuter_LE_StripNbr>>h11(1250,0,1250)","","");
c1->cd(12);
tree->Draw("fOuter_TE_StripNbr>>h12(1250,0,1250)","","");
c1->cd(13);
tree->Draw("fOuter_TE_StripNbr:fDetected_Position_Z>>h13","","colz");
c1->cd(14);
tree->Draw("fOuter_LE_StripNbr:fDetected_Position_X>>h14","","colz");
c1->cd(15);
tree->Draw("fRC_Theta[1]>>h15","","");
tree->Draw("fRC_Theta[1]>>h16","Strasse.GetOuterMultTEnergy()+Strasse.GetInnerMultTEnergy()","same");
h16->SetLineColor(2);
h16->Scale(1./4);;
c1->cd(16);
tree->Draw("fRC_Theta[0]>>h17","","");
tree->Draw("fRC_Theta[0]>>h18","Strasse.GetOuterMultTEnergy()+Strasse.GetInnerMultTEnergy()","same");
h18->SetLineColor(2);
h18->Scale(1./4);;


}
