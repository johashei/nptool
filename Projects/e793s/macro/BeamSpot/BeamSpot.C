/***************************************************************************
 *   Charlie Paxman (cp00474@surrey.ac.uk)                                 *
 ***************************************************************************/

//-------------------------------
//C++
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
//Root
#include <TVector3.h>
//NPTool
#include "NPEnergyLoss.h"
#include "NPReaction.h"
#include "NPSystemOfUnits.h"
using namespace std;
//-------------------------------

// TO DO:
// - Implement ROOT minimiser?
// - Too much repeptitive code! clean up, chuck in functions
// - Push project?

void BeamSpot(){
  vector <double> Xp, Yp, Zp, Ep;  //XYZE of detected particle. Constant.
  vector <int>    DetNum;          //Telescope number of hit
  vector <double> thetaLab, phiLab;
  double ELab = 0.0, Ex = 0.0;
  vector <double> Xd, Yd, Zd;      //Vector of particle direction. Calculated as Xp-Xb, Yp-Yb...
  ifstream MugastDataFile;

  /*** ITERATIVE GRID CONTROLS ***/
  /***** pos varied as offset ****/
  /**/ double xmin = +0.000;   /**/
  /**/ double xmax = +0.100;   /**/
  /**/ unsigned int xdiv = 10; /**/
  /**/                         /**/
  /**/ double ymin = -0.100;   /**/
  /**/ double ymax = +0.100;   /**/
  /**/ unsigned int ydiv = 10; /**/
  /**/                         /**/
  /**/ double zmin = -0.000;   /**/
  /**/ double zmax = +0.000;   /**/
  /**/ unsigned int zdiv =  1; /**/
  /**/                         /**/
  /***** thick varied as %ge *****/
  /**/ unsigned int tmin = 7;  /**/
  /**/ unsigned int tmax = 13; /**/
  /**/ double tmult = 0.1;     /**/
  /*******************************/

  // Calculate size of iteratve steps
  double xstp = (xmax-xmin)/ ((double) xdiv);
  double ystp = (ymax-ymin)/ ((double) ydiv);
  double zstp = (zmax-zmin)/ ((double) zdiv);
  cout << "Xstp = " << xstp << "   Ystp = " << ystp << "   Zstp = " << zstp << endl;

  // Vectors of the normal for each detector. UPDATE WITH NEW POSITIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  TVector3 MugastNormal1{ -0.453915, +0.455463, -0.765842};
  TVector3 MugastNormal2{ -0.642828, +0.000000, -0.766010};
  TVector3 MugastNormal3{ -0.454594, -0.450670, -0.768271};
  TVector3 MugastNormal4{ -0.002437, -0.638751, -0.769409};
  TVector3 MugastNormal5{ +0.452429, -0.454575, -0.767248};
  TVector3 MugastNormal7{ +0.443072, +0.443265, -0.779232};

  // Read in NPTool values
  NPL::EnergyLoss LightAl;
  NPL::EnergyLoss LightTarget;
  LightAl = NPL::EnergyLoss("proton_Al.G4table", "G4Table", 100);
  LightTarget = NPL::EnergyLoss("proton_CD2.G4table", "G4Table", 100);
  NPL::Reaction reaction;
  reaction.ReadConfigurationFile("../../Reaction/47Kdp_0keV.reaction");

  // Open and read event data file
  //MugastDataFile.open("XYZE_gammaGated_Full.txt", ios::in);
  MugastDataFile.open("XYZE_gammaGated_Run63.txt", ios::in);
  //MugastDataFile.open("XYZE_writeRun62_May04.txt", ios::in);
  //MugastDataFile.open("XYZE_writeRun63_May05.txt", ios::in);
  if(!MugastDataFile){
    cout << "ERROR: File not opened." << endl;
  }
  else{
    int numTemp = 0;
    double Xptemp = 0.0, Yptemp = 0.0, Zptemp = 0.0, Eptemp = 0.0; 
    while(true){ // Exit through break command
      MugastDataFile >> numTemp >> Xptemp >> Yptemp >> Zptemp >> Eptemp;
      DetNum.push_back(numTemp);
      Xp.push_back(Xptemp);
      Yp.push_back(Yptemp);
      Zp.push_back(Zptemp);
      Ep.push_back(Eptemp);
      
      Xptemp = Yptemp = Zptemp = Eptemp = 0.0;
      numTemp = 0;

      if(MugastDataFile.eof()) {break;}
    }
  }

  // Initilize various values and histograms
  int size = Xp.size();
  double tempTheta;
  TH1F *tempHist = new TH1F("tempHist","All Detectors", 40,-1.0,1.0);
  TH1F *tempMG1 = new TH1F("tempMG1","Individual MG#", 40,-1.0,1.0);
  TH1F *tempMG2 = new TH1F("tempMG2","Individual MG#", 40,-1.0,1.0);
  TH1F *tempMG3 = new TH1F("tempMG3","Individual MG#", 40,-1.0,1.0);
  TH1F *tempMG4 = new TH1F("tempMG4","Individual MG#", 40,-1.0,1.0);
  TH1F *tempMG5 = new TH1F("tempMG5","Individual MG#", 40,-1.0,1.0);
  TH1F *tempMG7 = new TH1F("tempMG7","Individual MG#", 40,-1.0,1.0);

  // ------------ ITERATE BEAM POSITION ------------ //
  for(int thickiter=tmin; thickiter<tmax+1; thickiter++){
    double TargetThickness = (thickiter * tmult) * 0.00476; //multiplier * original thickness
    //cout << " Thickness @ " << TargetThickness << endl;
    for(int xiter=0; xiter<xdiv+1; xiter++){
      //cout << " X @ " << xmin+(xiter*xstp) << endl;
      for(int yiter=0; yiter<ydiv+1; yiter++){
        //cout << " Y @ " << ymin+(yiter*ystp) << endl;
        for(int ziter=0; ziter<zdiv+1; ziter++){
          //cout << " Z @ " << zmin+(ziter*zstp) << endl;

  	  // Build beam spot vector for this iteration 
          TVector3 beamSpot{ xmin+(xiter*xstp), ymin+(yiter*ystp), zmin+(ziter*zstp) }; 

	  // Clean histograms
	  tempHist->Reset();
	  tempMG1->Reset();
	  tempMG2->Reset();
	  tempMG3->Reset();
	  tempMG4->Reset();
	  tempMG5->Reset();
	  tempMG7->Reset();

          // -------- ITERATE PARTICLE NUM -------- //
          for(int i=0; i<size-1; i++){
            // Build particle direction vector
	    TVector3 particleDir{ Xp[i] - beamSpot.X(),   //Xd
  	 	                  Yp[i] - beamSpot.Y(),   //Yd
			          Zp[i] - beamSpot.Z() }; //Zd
	  
	    tempTheta = ELab = Ex = 0.0;

	    switch(DetNum[i]){
	      case 1:
                tempTheta = particleDir.Angle(MugastNormal1);
	        break;
	      case 2:
                tempTheta = particleDir.Angle(MugastNormal2);
	        break;
	      case 3:
                tempTheta = particleDir.Angle(MugastNormal3);
	        break;
  	      case 4:
                tempTheta = particleDir.Angle(MugastNormal4);
	        break;
	      case 5:
                tempTheta = particleDir.Angle(MugastNormal5);
	        break;
	      case 7:
                tempTheta = particleDir.Angle(MugastNormal7);
	        break;
	      default:
	        cout << "ERROR! Invalid DetNum: " << DetNum[i] << " @" << i << endl;
	        return; // Exit code
	    }

	    //micrometer defined in NPSystemOfUnits.h
	    ELab = LightAl.EvaluateInitialEnergy(Ep[i], 0.4*micrometer, tempTheta); 
	    ELab = LightTarget.EvaluateInitialEnergy(ELab, 0.5*TargetThickness, 0.);

            // Change beam spot vector to beam direction vector
	    beamSpot.SetZ(1.0);
            Ex = reaction.ReconstructRelativistic( ELab, particleDir.Angle(beamSpot) );

	    // Fill Ex histograms
	    tempHist->Fill(Ex);
	    switch(DetNum[i]){
	      case 1:
                tempMG1->Fill(Ex); 
	        break;
	      case 2:
                tempMG2->Fill(Ex); 
	        break;
	      case 3:
                tempMG3->Fill(Ex); 
	        break;
  	      case 4:
                tempMG4->Fill(Ex); 
	        break;
	      case 5:
                tempMG5->Fill(Ex); 
	        break;
	      case 7:
                tempMG7->Fill(Ex); 
	        break;
	      default:
	        cout << "ERROR! Invalid DetNum: " << DetNum[i] << " @" << i << endl;
	        return; // Exit code
	    }
          }
	  // ------ END OF PARTICLE ITERATION ----- //
          // -------------------------------------- // 

	  // Initilise gaussian variable arrays
          double mean[7];
          double meanErr[7];
          double sigma[7];
          double sigmaErr[7];

	  // Draw Ex histograms
          TCanvas *c1 = new TCanvas("c1","Ex Histograms",20,20,1600,800);
          gStyle->SetOptStat(0);
          c1->Divide(2,1,0.005,0.005,0);
          c1->cd(1)->SetLeftMargin(0.15);
          c1->cd(1)->SetBottomMargin(0.15);
          gPad->SetTickx();
          gPad->SetTicky();
          c1->cd(2)->SetLeftMargin(0.15);
          c1->cd(2)->SetBottomMargin(0.15);
          gPad->SetTickx();
          gPad->SetTicky();
      
          c1->cd(1);
          tempMG1->SetMaximum(75.);
          tempMG1->GetXaxis()->SetTitle("Ex [MeV]");
          tempMG1->GetYaxis()->SetTitle("Counts");
      
          // ----- MG1 -----
          tempMG1->SetLineColor(kRed);
          tempMG1->SetFillStyle(3244);
          tempMG1->SetFillColor(kRed);
          tempMG1->Draw();
      
          tempMG1->Fit("gaus","WQ"); //add N to stop it drawing
          TF1 *gaus = (TF1*)tempMG1->GetListOfFunctions()->FindObject("gaus");
          mean[1] = gaus->GetParameter(1);
          meanErr[1] = gaus->GetParError(1);
          sigma[1] = gaus->GetParameter(2);
          sigmaErr[1] = gaus->GetParError(2);
          // ---------------
          // ----- MG2 -----
          tempMG2->SetLineColor(kOrange);
          tempMG2->SetFillStyle(3244);
          tempMG2->SetFillColor(kOrange);
          tempMG2->Draw("same");

          tempMG2->Fit("gaus","WQ"); //add N to stop it drawing
          gaus = (TF1*)tempMG2->GetListOfFunctions()->FindObject("gaus");
          mean[2] = gaus->GetParameter(1);
          meanErr[2] = gaus->GetParError(1);
          sigma[2] = gaus->GetParameter(2);
          sigmaErr[2] = gaus->GetParError(2);
          // ---------------
          // ----- MG3 -----
          tempMG3->SetLineColor(kGreen);
          tempMG3->SetFillStyle(3344);
          tempMG3->SetFillColor(kGreen);
          tempMG3->Draw("same");

          tempMG3->Fit("gaus","WQ"); //add N to stop it drawing
          gaus = (TF1*)tempMG3->GetListOfFunctions()->FindObject("gaus");
          mean[3] = gaus->GetParameter(1);
          meanErr[3] = gaus->GetParError(1);
          sigma[3] = gaus->GetParameter(2);
          sigmaErr[3] = gaus->GetParError(2);
          // ---------------
          // ----- MG4 -----
          tempMG4->SetLineColor(kTeal);
          tempMG4->SetFillStyle(3444);
          tempMG4->SetFillColor(kTeal);
          tempMG4->Draw("same");

	  tempMG4->Fit("gaus","WQ"); //add N to stop it drawing
          gaus = (TF1*)tempMG4->GetListOfFunctions()->FindObject("gaus");
          mean[4] = gaus->GetParameter(1);
          meanErr[4] = gaus->GetParError(1);
          sigma[4] = gaus->GetParameter(2);
          sigmaErr[4] = gaus->GetParError(2);
          // ---------------
          // ----- MG5 -----
	  tempMG5->SetLineColor(kBlue);
          tempMG5->SetFillStyle(3544);
          tempMG5->SetFillColor(kBlue);
          tempMG5->Draw("same");

          tempMG5->Fit("gaus","WQ"); //add N to stop it drawing
          gaus = (TF1*)tempMG5->GetListOfFunctions()->FindObject("gaus");
          mean[5] = gaus->GetParameter(1);
          meanErr[5] = gaus->GetParError(1);
          sigma[5] = gaus->GetParameter(2);
          sigmaErr[5] = gaus->GetParError(2);
          // ---------------
	  // ----- MG7 -----
	  tempMG7->SetLineColor(kViolet);
          tempMG7->SetFillStyle(3644);
          tempMG7->SetFillColor(kViolet);
	  tempMG7->Draw("same");

	  tempMG7->Fit("gaus","WQ"); //add N to stop it drawing
	  gaus = (TF1*)tempMG7->GetListOfFunctions()->FindObject("gaus");
          mean[6] = gaus->GetParameter(1);
          meanErr[6] = gaus->GetParError(1);
          sigma[6] = gaus->GetParameter(2);
          sigmaErr[6] = gaus->GetParError(2);
          // ---------------
          // Format legend
	  auto legend = new TLegend(0.15,0.7,0.35,0.9);
          legend->AddEntry(tempMG1,"MUGAST 1","f");
          legend->AddEntry(tempMG2,"MUGAST 2","f");
          legend->AddEntry(tempMG3,"MUGAST 3","f");
          legend->AddEntry(tempMG4,"MUGAST 4","f");
          legend->AddEntry(tempMG5,"MUGAST 5","f");
          legend->AddEntry(tempMG7,"MUGAST 7","f");
          legend->Draw();
          // ----- ALL -----
	  c1->cd(2);
          tempHist->GetXaxis()->SetTitle("Ex [MeV]");
          tempHist->GetYaxis()->SetTitle("Counts");
          tempHist->Draw();
          tempHist->Fit("gaus", "WQ");

	  gaus = (TF1*)tempHist->GetListOfFunctions()->FindObject("gaus");
          mean[0] = gaus->GetParameter(1);
          meanErr[0] = gaus->GetParError(1);
          sigma[0] = gaus->GetParameter(2);
          sigmaErr[0] = gaus->GetParError(2);
          // ---------------

	  // Caluclate metric
          double metric[7];
          for (int i=0; i<7; i++){
            metric[i] = pow(mean[i]-0.143,2) + pow(0.5 * sigma[i],2);
          }

	  // Write values to screen
          //cout << TargetThickness << "\t"
	  //     << xmin+(xiter*xstp) << "\t"
	  //     << ymin+(yiter*ystp) << "\t"
	  //     << zmin+(ziter*zstp) << "\t\t"
	  //     << metric[0] << "\t"
          //     << metric[1] << "\t"
          //     << metric[2] << "\t"
          //     << metric[3] << "\t"
          //     << metric[4] << "\t"
          //     << metric[5] << "\t"
          //     << metric[6] << endl;

	  cout << TargetThickness << "\t"
	       << xmin+(xiter*xstp) << "\t"
	       << ymin+(yiter*ystp) << "\t"
	       << zmin+(ziter*zstp) << "\t\t"
	       << mean[0] << "\t" << sigma[0] << "\t"
	       << mean[1] << "\t" << sigma[1] << "\t"
	       << mean[2] << "\t" << sigma[2] << "\t"
	       << mean[3] << "\t" << sigma[3] << "\t"
	       << mean[4] << "\t" << sigma[4] << "\t"
	       << mean[5] << "\t" << sigma[5] << "\t"
	       << mean[6] << "\t" << sigma[6] << "\t"
               << endl;


	  // Write histograms to file
          string fileOut = "./Histograms/targetIter";
          fileOut += to_string(TargetThickness);
          fileOut += "_x";
          fileOut += to_string(xmin+(xiter*xstp));
          fileOut += "_y";
          fileOut += to_string(ymin+(yiter*ystp));
          fileOut += "_z";
          fileOut += to_string(zmin+(ziter*zstp));
          fileOut += ".pdf";

          c1->SaveAs(fileOut.c_str());
        } //iterate Z
      } //iterate Y
    } //iterate X
  } //iterate thickness
  //cout << " --- COMPLETE --- " << endl;
  // ----------------------------------------------- //
}
