// Contain global variable declaration, comment and option
     //cout << "fELabFront " << fELabFront << " fELabBack " << fELabBack << " fELabHalf " << fELabHalf << " half calcul " << (fELabFront+fELabBack)/2 << endl;
#include "EnergyCalibrator.h"

void AutoCalibration(int Telescope_Start, int Telescope_End, Double_t Dispersion_Limit, Double_t Al_InitialThickness)
{
	for(int i = Telescope_Start ; i<=Telescope_End ;i++)
	{
		AlThickness = Al_InitialThickness*micrometer ; // initial thickness that will be adjusted
    cout << "Al thickness ini = " << AlThickness << endl;
		SiThickness = 0.0*micrometer ; // not taken into account
		double Al_step = 0.01*micrometer;
		int step_limit = 10; 
		int k = 0 ; 

		Telescope_Number = i ;
		// Create a folder to Hold all the file from calibration
		ostringstream FolderName;
		FolderName << Experiment << "_" << xy << "_D" << Telescope_Number << "_E";
		main_name = FolderName.str() ;
		TString make_folder = "mkdir -p ./Calibration207Bi/" + main_name ;   
		folder = "./Calibration207Bi/" + FolderName.str() ;
		int returnValue1 = system(make_folder);
		int returnValue2 = system(make_folder+"/peaks");
		int returnValue3 = system(make_folder+"/dispersion");
		int returnValue4 = system(make_folder+"/latex");
		int returnValue5 = system(make_folder+"/latex/pictures");

		// open the ROOT file alpha to process
		TString path_histo  = "./Histograms/";
		TString inFileName = frun_207Bi;
		inFile = new TFile(path_histo + inFileName +"_RawDSSSDHistos.root");

    // If pedestal is obtained with pulser run
    if(Pedestals_pulser) {
      // Open ROOT file containing calibration coefficients obtained for the pulser runs
      TString path_calib  = "./calibs/";
      TString inFileName_pf = frun_pulserFront;
      TString inFileName_pb = frun_pulserBack;
      f_pf = new TFile(path_calib + inFileName_pf + "_pedestals_front.root");
      f_pb = new TFile(path_calib + inFileName_pb + "_pedestals_back.root");
    }

    // Run EnergyCalibrator
		EnergyCalibrator();
		bool check1=false,check2=false;
		//while( !(mean_extrapolation <0.1 && mean_extrapolation >-0.1 ) && k < step_limit )
		while( !(mean_extrapolation <Dispersion_Limit & mean_extrapolation >-Dispersion_Limit) && k < step_limit ) // Dispersion_Limit equivalent to sigma
		{
      if(mean_extrapolation < 0)
			{
//				if(xy=="FRONT") {
					AlThickness -= Al_step;
          cout << "extrapo < 0 Al thickness step- = " << AlThickness << endl;
/*        }
				else if(xy=="BACK") {
					AlThickness += Al_step;
          cout << "extrapo > 0 Al thickness step+ = " << AlThickness << endl;
        }*/
				check1=true;
			}

			else if (mean_extrapolation > 0)
			{
//				if(xy=="FRONT") {
					AlThickness += Al_step;
          cout << " extrapo > 0 Al thickness step+ = " << AlThickness << endl;
/*        }
				else if(xy=="BACK") {
					AlThickness -= Al_step;
          cout << " extrapo < 0 Al thickness step+ = " << AlThickness << endl;
        }*/
				check2=true;
			}

			if(check1&&check2)
			{
				Al_step=Al_step/10.;
				check1=false;check2=false;
        cout << "check1&&check2" << endl;
			}
			latex_file.close();
			EnergyCalibrator();

			cout << " Iteration Results: Al Thickness: " << AlThickness/micrometer << "um | Mean Extrapolation  "  << mean_extrapolation << "Chan. "<< endl ;
      cout << endl;

			k++;
      cout << "k = " << k << endl;
		}

		LatexSummaryEnder();

		delete Buffer;
		delete Source_branching_ratio;
		delete Source_E;
		delete Source_Sig;
		delete energyFront;
		delete errorsFront;
		delete energyBack;
		delete errorsBack;

	}

	return;
}


/////////////////////////////
void DefineSource(TString sourceName, Double_t angle)
//void DefineSource(TString sourceName)
{
	if(sourceName=="3 alphas")
	{
		NumberOfIsotope = 3 ;
		energyFront = new Double_t[NumberOfIsotope]; errorsFront = new Double_t[NumberOfIsotope];
		energyBack = new Double_t[NumberOfIsotope]; errorsBack = new Double_t[NumberOfIsotope];

		/// Information used in the summary
		Source_Number_Peak = 8;
		Source_isotope = new TString[Source_Number_Peak] ;Source_E = new Double_t[Source_Number_Peak] ; Source_Sig = new Double_t[Source_Number_Peak] ; Source_branching_ratio = new Double_t[Source_Number_Peak] ;

		// 239Pu
		Source_isotope[0]="$^{239}$Pu"; Source_E[0]   = 5.15659 ; Source_Sig[0] = 0.00014 ; Source_branching_ratio[0] = 70.77 ;
		Source_isotope[1]="$^{239}$Pu"; Source_E[1]   = 5.1443 ; Source_Sig[1] = 0.00014 ; Source_branching_ratio[1] = 17.11 ;
		Source_isotope[2]="$^{239}$Pu"; Source_E[2]   = 5.1055  ; Source_Sig[2] = 0.00014 ; Source_branching_ratio[2] = 11.94 ;

		// 241Am
		Source_isotope[3]="$^{241}$Am"; Source_E[3]   = 5.48556 ; Source_Sig[3] = 0.00012 ; Source_branching_ratio[3] = 84.8 ;
		Source_isotope[4]="$^{241}$Am"; Source_E[4]   = 5.44280 ; Source_Sig[4] = 0.00012 ; Source_branching_ratio[4] = 13.1 ;
		Source_isotope[5]="$^{241}$Am"; Source_E[5]   = 5.388   ; Source_Sig[5] = 0.00012 ; Source_branching_ratio[5] = 1.66 ;

		// 244Cm
		Source_isotope[6]="$^{244}$Cm"; Source_E[6]   = 5.80477 ; Source_Sig[6] = 0.00005 ; Source_branching_ratio[6] = 76.90 ;
		Source_isotope[7]="$^{244}$Cm"; Source_E[7]   = 5.76264 ; Source_Sig[7] = 0.00005 ; Source_branching_ratio[7] = 23.10 ;

		// Corrected value of main peak used in the fit to have a reasonable chi2
		Double_t sig_value = 0.001;
		//Double_t sig_value = 0.1;
		Double_t alpha1_Sig = sig_value ; Double_t alpha2_Sig = sig_value ; Double_t alpha3_Sig = sig_value ; 

    //cout << "Etheo1 ini = " <<  Source_E[0] << endl; 
		
    Double_t alpha1_E , alpha2_E , alpha3_E;

		alpha1_E = EL_Al.Slow(	Source_E[0]*MeV , // Energy of the detected particle
				AlThickness	    , // Target Thickness at 0 degree
				angle			          ) ; // particle angle
    //cout << "particle angle = " << angle << endl;
    //cout << "Etheo1 Al = " << alpha1_E << endl;

		alpha1_E = EL_Si.Slow(	alpha1_E*MeV    , // Energy of the detected particle
				SiThickness	    , // Target Thickness at 0 degree
				0			          ) ; // particle angle
    //cout << "Etheo1 Si = " << alpha1_E << endl;

		alpha2_E = EL_Al.Slow(	Source_E[3]*MeV , // Energy of the detected particle
				AlThickness	    , // Target Thickness at 0 degree
				angle			          ) ; // particle angle

		alpha2_E = EL_Si.Slow(	alpha2_E*MeV    , // Energy of the detected particle
				SiThickness	    , // Target Thickness at 0 degree
				0     		      ) ; // particle angle

		alpha3_E = EL_Al.Slow(	Source_E[6]*MeV , // Energy of the detected particle
				AlThickness	    , // Target Thickness at 0 degree
				angle			          ) ; // particle angle

		alpha3_E = EL_Si.Slow(	alpha3_E*MeV    , // Energy of the detected particle
				SiThickness   	, // Target Thickness at 0 degree
				0     		      ) ; // particle angle           


		// front and back are in a same order        
		energyFront[0] = alpha1_E   ; energyFront[1] = alpha2_E   ; energyFront[2] = alpha3_E   ;
		errorsFront[0] = alpha1_Sig ; errorsFront[1] = alpha2_Sig ; errorsFront[2] = alpha3_Sig ;

    if (Pedestals_pulser) {
      energyBack[0] = alpha1_E   ; energyBack[1] = alpha2_E   ; energyBack[2] = alpha3_E   ;
      errorsBack[0] = alpha1_Sig ; errorsBack[1] = alpha2_Sig ; errorsBack[2] = alpha3_Sig ;
    }
    else {
      energyBack[0] = alpha3_E   ; energyBack[1] = alpha2_E   ; energyBack[2] = alpha1_E   ;
      errorsBack[0] = alpha3_Sig ; errorsBack[1] = alpha2_Sig ; errorsBack[2] = alpha1_Sig ;
		}
	}

  else if(sourceName=="207Bi")
	{
		NumberOfPeaks = 6;
		energyFront = new Double_t[NumberOfPeaks]; errorsFront = new Double_t[NumberOfPeaks];
		energyBack = new Double_t[NumberOfPeaks]; errorsBack = new Double_t[NumberOfPeaks];

		/// Information used in the summary
		Source_Number_Peak = 6;
		Source_isotope = new TString[Source_Number_Peak] ;Source_E = new Double_t[Source_Number_Peak] ; Source_Sig = new Double_t[Source_Number_Peak] ; Source_branching_ratio = new Double_t[Source_Number_Peak] ;

    // E in MeV
		// e-(CE K) 481.6935 keV
		Source_isotope[0]="e(CEK)482"; Source_E[0]   = 481.6935e-3 ; Source_Sig[0] = 0.0021e-3 ; Source_branching_ratio[0] = 1.537 ;
    // e-(CE L) 553.8372 keV
    Source_isotope[1]="e(CEL)554"; Source_E[1]   = 553.8372e-3 ; Source_Sig[1] = 0.0021e-3 ; Source_branching_ratio[1] = 0.442 ;
    // e-(CE M) 565.8473 keV
    Source_isotope[2]="e(CEM)566"; Source_E[2]   = 565.8473e-3  ; Source_Sig[2] = 0.0021e-3 ; Source_branching_ratio[2] = 0.111 ;
		// e-(CE K) 975.651 keV
		Source_isotope[3]="e(CEK)976"; Source_E[3]   = 975.651e-3 ; Source_Sig[3] = 0.003e-3 ; Source_branching_ratio[3] = 7.08 ;
    // e-(CE L) 1047.795 keV
    Source_isotope[4]="e(CEL)1048"; Source_E[4]   = 1047.795e-3 ; Source_Sig[4] = 0.003e-3 ; Source_branching_ratio[4] = 1.84 ;
    // e-(CE M) 1059.805
    Source_isotope[5]="$e(CEM)1060"; Source_E[5]   = 1059.805e-3   ; Source_Sig[5] = 0.003e-3 ; Source_branching_ratio[5] = 0.44 ;

		// Corrected value of main peak used in the fit to have a reasonable chi2
		Double_t sig_value = 0.001;
		//Double_t sig_value = 0.1;
		Double_t peak1_Sig = sig_value ; Double_t peak2_Sig = sig_value ; Double_t peak3_Sig = sig_value ; 
		Double_t peak4_Sig = sig_value ; Double_t peak5_Sig = sig_value ; Double_t peak6_Sig = sig_value ; 

    //cout << "Etheo1 ini = " <<  Source_E[0] << endl; 
		
    Double_t peak1_E , peak2_E , peak3_E;
    Double_t peak4_E , peak5_E , peak6_E;

/*		alpha1_E = EL_Al.Slow(	Source_E[0]*MeV , // Energy of the detected particle
				AlThickness	    , // Target Thickness at 0 degree
				angle			          ) ; // particle angle
    //cout << "particle angle = " << angle << endl;
    //cout << "Etheo1 Al = " << alpha1_E << endl;

		alpha1_E = EL_Si.Slow(	alpha1_E*MeV    , // Energy of the detected particle
				SiThickness	    , // Target Thickness at 0 degree
				0			          ) ; // particle angle
    //cout << "Etheo1 Si = " << alpha1_E << endl;

		alpha2_E = EL_Al.Slow(	Source_E[3]*MeV , // Energy of the detected particle
				AlThickness	    , // Target Thickness at 0 degree
				angle			          ) ; // particle angle

		alpha2_E = EL_Si.Slow(	alpha2_E*MeV    , // Energy of the detected particle
				SiThickness	    , // Target Thickness at 0 degree
				0     		      ) ; // particle angle

		alpha3_E = EL_Al.Slow(	Source_E[6]*MeV , // Energy of the detected particle
				AlThickness	    , // Target Thickness at 0 degree
				angle			          ) ; // particle angle

		alpha3_E = EL_Si.Slow(	alpha3_E*MeV    , // Energy of the detected particle
				SiThickness   	, // Target Thickness at 0 degree
				0     		      ) ; // particle angle           
*/

		// front and back are in a same order        
		energyFront[0] = peak1_E   ; energyFront[1] = peak2_E   ; energyFront[2] = peak3_E   ;
		energyFront[3] = peak4_E   ; energyFront[4] = peak5_E   ; energyFront[5] = peak6_E   ;
		errorsFront[0] = peak1_Sig ; errorsFront[1] = peak2_Sig ; errorsFront[2] = peak3_Sig ;
		errorsFront[3] = peak4_Sig ; errorsFront[4] = peak5_Sig ; errorsFront[5] = peak6_Sig ;

    energyBack[0] = peak1_E   ; energyBack[1] = peak2_E   ; energyBack[2] = peak3_E   ;
		energyBack[3] = peak4_E   ; energyBack[4] = peak5_E   ; energyBack[5] = peak6_E   ;
		errorsBack[0] = peak1_Sig ; errorsBack[1] = peak2_Sig ; errorsBack[2] = peak3_Sig ;
		errorsBack[3] = peak4_Sig ; errorsBack[4] = peak5_Sig ; errorsBack[5] = peak6_Sig ;


	}

	return;  
}


/////////////////////////////
void EnergyCalibrator()
{
	// Set-up the root Style
	gStyle->SetOptTitle();
	gStyle->SetOptTitle();

	//DefineSource();

	TString str;
	TString str1;
	TString str2;
	TString strbuff;
	TString strbuff2;
	TString fname;
	TString fname2;
	TString fname3;
	TString fname4;
	TString hname;

	LatexSummaryHeader(xy);

	// Clear everything
	BadStrip.clear() ;
	sigma_fit = new TH1F("Sigma", "Sigma from fit (channel)", 80, 0,10);
	Dispersion= new TH1F("Dispersion", "Dispersion from Zero Extrapolation (channel)", 40, -20,20);
	ZeroDispersion = new TGraph(); // 16 ?
	coeff_a = new TGraph(); // 16 ?
	coeff_b = new TGraph(); // 16 ?

	ostringstream number ;
	number << Telescope_Number  ;
	CurrentTelescope = Telescope_Number ;
	if (xy == "FRONT"){ 
		//////// Input Files ///////////
		str = "hCOMPTONTELESCOPE_D"+number.str()+"_FRONT_E";
		/////// Output Files ///////////
		str1 = "Cal_Str_FRONT_E_COMPTONTELESCOPE_"+number.str();
	}	 
	else if (xy == "BACK"){ 	
		//////// Input Files ///////////
		str = "hCOMPTONTELESCOPE_D"+number.str()+"_BACK_E";
		/////// Output Files ///////////
		str1 = "Cal_Str_BACK_E_COMPTONTELESCOPE_"+number.str();
	}	 
	else {cout<< "ERROR FOR FRONT or BACK PARAMETER"<< endl;}

	fname =  folder + "/peaks/" + str1 + ".peak";
	peaks_file.open( ( (string)fname ).c_str() );

	fname2 = folder + "/" + str1 + ".cal";
	calib_file.open( ( (string)fname2 ).c_str() );

	fname3 = folder + "/" + str1 + ".dispersion";
	dispersion_file.open( ( (string)fname3 ).c_str() );

  fname4 = folder + "/" + str1 + ".calError";
  calibError_file.open( ( (string)fname4 ).c_str() );

	Tsummary = new TCanvas((TString)("Telescope"+number.str()+"Summary"), (TString)("Telescope "+number.str()+" Summary"), 700, 700);
	Tsummary->Divide(2,3);
	Buffer  = new TCanvas((TString)("Buffer"), (TString)("Buffer"), 10, 10);
	Buffer->cd(1);

	for (Int_t j = Strip_Start-1; j < Strip_End; j++)
	{
    
		CurrentStrip=j+1;
		number.seekp(0);
		number << j+1;

    // Get particle angle and determine the 3 alpha peaks energies after energy loss

    Double_t particle_angle = 0;

/*    if (xy == "FRONT") {
      // angle = 0
      Double_t particle_angle = 0;
      DefineSource(sourceName, particle_angle);
    }
    
    else if (xy == "BACK"){
      // angle = 0
      Double_t particle_angle = 0;
      DefineSource(sourceName, particle_angle);
    }
*/
    DefineSource(sourceName, particle_angle);

		///// Get the alpha histogram of det i and strip j /////
		hname = str+number.str();
		TH1F *hist = (TH1F*) inFile->Get(((string)hname).c_str());

		// Prevent rebinning in Pedestal
		TH1F *histSource = (TH1F*)hist->Clone();

    if (Pedestals_pulser) {
    ///// Get offset graph from pulser calibration
      //cout << "///////////// pedestal pulser /////////////" << endl;
      if (xy == "FRONT") {
        // get graph
        TGraph *gr_pof = (TGraph*) f_pf->Get("grOffset");
        // get total number of strips, absolute_strip and offset
        Int_t nstrip_pof = gr_pof->GetN();
        Double_t *x_pof = gr_pof->GetX();
        Double_t *y_pof = gr_pof->GetY();
        // associate current strip number to strip number in root file
        Int_t absolute_CurrentStrip = 2*(CurrentTelescope-1)*Strip_End + CurrentStrip;
        cout << "Telescope = " << CurrentTelescope << " strip = " << CurrentStrip << /*" absolute strip equiv = " << absolute_CurrentStrip <<*/ endl;
        // loop on pedestals offset and get offset for the current strip
        for (Int_t m = 0; m < nstrip_pof; m++) {
          if (x_pof[m] == absolute_CurrentStrip) {
            pedestal = y_pof[m];
//          cout << "strip equal ? " << x_pof[m] << " = " << absolute_CurrentStrip << " ?" << endl;
            //cout << "pedestal = " << pedestal << endl;
          }
        }
      }
      if (xy == "BACK") {
        // get graph
        TGraph *gr_pob = (TGraph*) f_pb->Get("grOffset");
        // get total number of strips, absolute_strip and offset
        Int_t nstrip_pob = gr_pob->GetN();
        Double_t *x_pob = gr_pob->GetX();
        Double_t *y_pob = gr_pob->GetY();
        // associate current strip number to strip number in root file
        Int_t absolute_CurrentStrip = (2*(CurrentTelescope-1)+1)*Strip_End + CurrentStrip;
        cout << "Telescope = " << CurrentTelescope << " strip = " << CurrentStrip /*<< " absolute strip equiv = " << absolute_CurrentStrip*/ << endl;
        // loop on pedestals offset and get offset for the current strip
        for (Int_t m = 0; m < nstrip_pob; m++) {
          if (x_pob[m] == absolute_CurrentStrip) {
            pedestal = y_pob[m];
//          cout << "strip equal ? " << x_pob[m] << " = " << absolute_CurrentStrip << " ? " << endl;
            //cout << "pedestal = " << pedestal << endl;
          }
        }
      }
    
      Alpha(histSource,
          xy,
          pedestal);
    }

    else {
      Alpha(histSource,
          xy,
          Pedestals(hist));
    }

	//		Tsummary->WaitPrimitive();
		histSource->GetXaxis()->SetRangeUser(0, 1024);
/*		if(xy == "FRONT") 	   histSource->GetXaxis()->SetRangeUser(0, 1024);
		//if(xy == "FRONT") 	   histSource->GetXaxis()->SetRangeUser(500, 3000);
		else if(xy == "BACK") histSource->GetXaxis()->SetRangeUser(0, 1024);
*/
/*		if(j == 127) 
		{ TH1F histSource67 = TH1F(*histSource);
			Tsummary->cd(2);
			histSource67.SetStats(true);
			histSource67.SetTitle("Raw Spectrum of strip 128 with gaussian fit");

			histSource67.Draw("A");
			Buffer->cd(1);
		}*/
	}

	Tsummary->cd(1); 
	mean_extrapolation = ZeroDispersion->GetMean(2);
  cout << "mean_extrapolation " << mean_extrapolation << endl;
	ZeroDispersion->SetMaximum(mean_extrapolation+10);ZeroDispersion->SetMinimum(mean_extrapolation-10);
	ZeroDispersion->SetTitle("Scattered plot of zero extrapolation dispersion : Ped.+b/a");
	ZeroDispersion->SetMarkerStyle(2);
	ZeroDispersion->Draw("ap");
	//Draw the mean line
	TLine mean_line = TLine(0, ZeroDispersion->GetMean(2), 17.5, ZeroDispersion->GetMean(2) );
	mean_line.Draw("");

	Tsummary->cd(3); 
	Dispersion->Draw();

	Tsummary->cd(4); 
	sigma_fit->Draw();

	Tsummary->cd(5) ;
	coeff_a->SetMarkerStyle(2);
	coeff_a->SetMaximum(coeff_a->GetMean(2)+0.0031);coeff_a->SetMinimum(coeff_a->GetMean(2)-0.0031);
	coeff_a->SetTitle("Gain a (MeV/channel)");
	coeff_a->Draw("ap");

	Tsummary->cd(6);
	coeff_b->SetMaximum(coeff_b->GetMean(2)+0.2);coeff_b->SetMinimum(coeff_b->GetMean(2)-0.2);
	coeff_b->SetMarkerStyle(2);
	coeff_b->SetTitle("Offset b (MeV)");
	coeff_b->Draw("ap");

	TString filename = Tsummary->GetName();
	Tsummary->SaveAs(filename+".pdf");
	Tsummary->Close();
	int returnValue6 = system("mv "+filename+".pdf ./" + folder + "/latex/pictures");

	peaks_file.close();
	calib_file.close();
	dispersion_file.close();
  calibError_file.close();

	LatexSummaryTelescope();
	delete Tsummary   ;
	delete sigma_fit  ;
	delete Dispersion ;

	Buffer->Close();   
}

/////////////////////////////////
Double_t Pedestals(TH1F *hist)
{

	if(Pedestals_Aligned)
		return 0 ;
		//return 8192 ;

	else
	{
    TF1 *gauss=new TF1("gauss","gaus",0,1024);

		hist->SetAxisRange(7800,8500);

		///// Peak search /////
		TSpectrum *s = new TSpectrum(2,1);
		Int_t nfound =0;
		nfound = s->Search(hist,2," ");

		Float_t *xpeaks = s->GetPositionX();

		Float_t linf =0, lsup =0; 
		Double_t sum=0, mean=0, sigma=0;

		if(nfound != 1 ) 
			cout << "########   PROBLEM Nfound != NAsked !  ########   " << hist->GetName() <<"  Nfound:"<<nfound<<endl;

		else {
			linf = xpeaks[0]-10;
			lsup = xpeaks[0]+10; 
			gauss=new TF1("gauss","gaus",linf,lsup); 
			gauss->SetRange(linf,lsup);
			hist->Fit(gauss,"RQ");

			sum = gauss->GetParameter(0);
			mean = gauss->GetParameter(1);
			sigma = gauss->GetParameter(2);

			if(sigma > 3)
				BadStrip[CurrentStrip] += " Alpha peak to large;" ;

		}

		delete s; delete gauss;
		return (mean) ;
	}
}

/////////////////////////////////
void Alpha(TH1F *hist, TString xy, Double_t Pedestal)
{

	hist->GetXaxis()->SetRangeUser(0,1024);
/*	if(xy == "FRONT") 		hist->GetXaxis()->SetRangeUser(500,4000);
	else if(xy == "BACK") 	hist->GetXaxis()->SetRangeUser(500,4000);
*/
	if(!Finder(hist, xy, mean, sigma )) cout << "On "<< hist->GetName() << endl ;

	// Fit
	TGraphErrors* gr_MM= new TGraphErrors(7); // 6 peaks + pedestal

	if(method == "ZeroForce")
	{
		a = Calib_ZeroForceMethod((string)xy ,gr_MM, Pedestal, mean, sigma);
		b = -Pedestal*a ;
	}

	else if(method == "ZeroExtrapolation")
	{
		Calib_ZeroExtrapolationMethod(hist, (string)xy, gr_MM, Pedestal, mean, error_mean, sigma, a, b);
	}

}

/////////////////////////////////
bool Finder(TH1F *h, TString xy, Double_t *mean, Double_t *sigma)
{

	/////////////////////////////////////////////////
	//						                                 //
	//	           PEAK  FINDER		                 //
	//						                                 //
	/////////////////////////////////////////////////

  // clone histo for later drawing
  TH1F *h2 = (TH1F*) h->Clone();

  if (sourceName == "3 alphas")
  {
    for(int k=0; k<3; k++)
    {
      mean[k]=0;
      error_mean[k]=0;
      sigma[k]=0;
    }

    Double_t resolsig=5;
    Float_t  resolsigTSpec=1;
    Double_t seuil=0.3;
    Int_t npeaks=8;   // maximum number of peak that can be found
    //Int_t npeaks=5;   // maximum number of peak that can be found

    //////// Peak finder

    TSpectrum *s = new TSpectrum(npeaks,resolsigTSpec);

    Int_t nfound = s->Search(h,resolsig,"new",seuil);
    Float_t *xpeaks = s->GetPositionX();

    /// Sort in growing order the array

    if(nfound>1)
    {
      for(Int_t p=0;p<nfound;p++)
      {
        for(Int_t i=0;i<nfound-1;i++)
        {
          if(xpeaks[i]>xpeaks[i+1])
          {
            Float_t varia=xpeaks[i];
            xpeaks[i]=xpeaks[i+1];
            xpeaks[i+1]=varia;
          }	  
        }
      }
    }

    Float_t linf=0, lsup=0; 

    // If 3 peak found
    if(nfound == 3)
    {
      for (Int_t p=0;p<nfound;p++)
      {
        if (xy == "FRONT") {
          if (CurrentTelescope == 2) {
            if (CurrentStrip == 1 ) {
              linf = xpeaks[p]-2;
              lsup = xpeaks[p]+6;
            }
            else {
              linf = xpeaks[p]-4;
              lsup = xpeaks[p]+8;
            }
          }
          else {
            linf = xpeaks[p]-4;
            lsup = xpeaks[p]+8;
          }
        }

        else if (xy == "BACK") {
          if (CurrentTelescope == 4) {
            if (CurrentStrip == 1 ) {
              linf = xpeaks[p]-13;
              lsup = xpeaks[p]+20;
            }
            else if (CurrentStrip == 16) {
              linf = xpeaks[p]-13;
              lsup = xpeaks[p]+20;
            }
            else {
              //linf = xpeaks[p]-2;
              linf = xpeaks[p]-7;
              //lsup = xpeaks[p]+10;
              lsup = xpeaks[p]+10;
            }
          }
          else {
            linf = xpeaks[p]-4;
            lsup = xpeaks[p]+8;
          }
        }

        /*			if(xy == "X")
                {			
                linf = xpeaks[p]-2;
                lsup = xpeaks[p]+8;


                else if (xy == "Y")
                {			
                linf = xpeaks[p]-8;
                lsup = xpeaks[p]+2;
                }
                */

        // fit
        TF1 *gauss = new TF1("gauss","gaus",linf,lsup);
        h->Fit(gauss,"RQ");
        mean[p] = gauss->GetParameter(1);
        error_mean[p] = gauss->GetParError(1);
        //cout << "mean and error " << mean[p] << "\t" << error_mean[p] << endl; 
        sigma[p]= gauss->GetParameter(2);

        sigma_fit->Fill(gauss->GetParameter(2));

        // save graph and fit in a root file
        TFile *outFile_FitHistoAlpha = new TFile("/projet/astronuc/31P3Het/nptool/Projects/nsi91/W1_Analysis/raw/graphs_FitAlpha.root", "update");
        h->SetNameTitle(Form("grFitAlpha_D%d_Strip0%d_peak%d", CurrentTelescope, CurrentStrip, p+1), Form("D%d_Strip0%d_peak%d", CurrentTelescope, CurrentStrip, p+1));
        if (xy == "FRONT" && CurrentTelescope == 2 && CurrentStrip == 1) {
          h->GetXaxis()->SetRangeUser(800, 1000);
        }
        else {
          h->GetXaxis()->SetRangeUser(1400, 2300);
        }
        h->Write();
        //outFile_FitHistoAlpha->Close();


      }
      }

      if(nfound!=3)
      {
        ostringstream numP;
        numP << nfound ;

        BadStrip[CurrentStrip] += " " + numP.str() + " peak(s) found;" ;

        for (Int_t p=0;p<3;p++)
        {
          cout << "attention, nombre de pics different de 3!!!" ;
          mean[p]=-1;
          sigma[p]=-1;
          return false ;
        }
      }

      // save raw histo
      TFile *outFile_RawHistoAlpha = new TFile("/projet/astronuc/31P3Het/nptool/Projects/nsi91/W1_Analysis/raw/graphs_rawAlpha.root", "update");
      h2->SetNameTitle(Form("grRawAlpha_D%d_Strip0%d", CurrentTelescope, CurrentStrip), Form("D%d_Strip0%d", CurrentTelescope, CurrentStrip));
      if (xy == "FRONT" && CurrentTelescope == 2 && CurrentStrip == 1) {
        h2->GetXaxis()->SetRangeUser(800, 1000);
      }
      else {
        h2->GetXaxis()->SetRangeUser(1400, 2300);
      }
      h2->Write();
      //outFile_FitHistoAlpha->Close();


      return true ;
    }

  else if (sourceName == "207Bi") 
  {
    for(int k=0; k<6; k++)
    {
      mean[k]=0;
      error_mean[k]=0;
      sigma[k]=0;
    }

    Double_t resolsig=5;
    Float_t  resolsigTSpec=1;
    Double_t seuil=0.3;
    Int_t npeaks=6;   // maximum number of peak that can be found

    //////// Peak finder

    TSpectrum *s = new TSpectrum(npeaks,resolsigTSpec);

    Int_t nfound = s->Search(h,resolsig,"new",seuil);
    Float_t *xpeaks = s->GetPositionX();

    /// Sort in growing order the array

    if(nfound>1)
    {
      for(Int_t p=0;p<nfound;p++)
      {
        for(Int_t i=0;i<nfound-1;i++)
        {
          if(xpeaks[i]>xpeaks[i+1])
          {
            Float_t varia=xpeaks[i];
            xpeaks[i]=xpeaks[i+1];
            xpeaks[i+1]=varia;
          }	  
        }
      }
    }

    Float_t linf=0, lsup=0; 

    // If 3 peak found
    if(nfound == 4)
    {
      for (Int_t p=0;p<nfound;p++)
      {
        if (xy == "FRONT") {
          if (CurrentTelescope == 2) {
            if (CurrentStrip == 1 ) {
              linf = xpeaks[p]-2;
              lsup = xpeaks[p]+6;
            }
            else {
              linf = xpeaks[p]-4;
              lsup = xpeaks[p]+8;
            }
          }
          else {
            linf = xpeaks[p]-4;
            lsup = xpeaks[p]+8;
          }
        }

        else if (xy == "BACK") {
          if (CurrentTelescope == 4) {
            if (CurrentStrip == 1 ) {
              linf = xpeaks[p]-13;
              lsup = xpeaks[p]+20;
            }
            else if (CurrentStrip == 16) {
              linf = xpeaks[p]-13;
              lsup = xpeaks[p]+20;
            }
            else {
              //linf = xpeaks[p]-2;
              linf = xpeaks[p]-7;
              //lsup = xpeaks[p]+10;
              lsup = xpeaks[p]+10;
            }
          }
          else {
            linf = xpeaks[p]-4;
            lsup = xpeaks[p]+8;
          }
        }

        /*			if(xy == "X")
                {			
                linf = xpeaks[p]-2;
                lsup = xpeaks[p]+8;


                else if (xy == "Y")
                {			
                linf = xpeaks[p]-8;
                lsup = xpeaks[p]+2;
                }
                */

        // fit
        TF1 *gauss = new TF1("gauss","gaus",linf,lsup);
        h->Fit(gauss,"RQ");
        mean[p] = gauss->GetParameter(1);
        error_mean[p] = gauss->GetParError(1);
        //cout << "mean and error " << mean[p] << "\t" << error_mean[p] << endl; 
        sigma[p]= gauss->GetParameter(2);

        sigma_fit->Fill(gauss->GetParameter(2));

        // save graph and fit in a root file
        TFile *outFile_FitHistoAlpha = new TFile("/projet/astronuc/31P3Het/nptool/Projects/nsi91/W1_Analysis/raw/graphs_FitAlpha.root", "update");
        h->SetNameTitle(Form("grFitAlpha_D%d_Strip0%d_peak%d", CurrentTelescope, CurrentStrip, p+1), Form("D%d_Strip0%d_peak%d", CurrentTelescope, CurrentStrip, p+1));
        if (xy == "FRONT" && CurrentTelescope == 2 && CurrentStrip == 1) {
          h->GetXaxis()->SetRangeUser(800, 1000);
        }
        else {
          h->GetXaxis()->SetRangeUser(1400, 2300);
        }
        h->Write();
        //outFile_FitHistoAlpha->Close();


      }
      }

      if(nfound!=3)
      {
        ostringstream numP;
        numP << nfound ;

        BadStrip[CurrentStrip] += " " + numP.str() + " peak(s) found;" ;

        for (Int_t p=0;p<3;p++)
        {
          cout << "attention, nombre de pics different de 3!!!" ;
          mean[p]=-1;
          sigma[p]=-1;
          return false ;
        }
      }

      // save raw histo
      TFile *outFile_RawHistoAlpha = new TFile("/projet/astronuc/31P3Het/nptool/Projects/nsi91/W1_Analysis/raw/graphs_rawAlpha.root", "update");
      h2->SetNameTitle(Form("grRawAlpha_D%d_Strip0%d", CurrentTelescope, CurrentStrip), Form("D%d_Strip0%d", CurrentTelescope, CurrentStrip));
      if (xy == "FRONT" && CurrentTelescope == 2 && CurrentStrip == 1) {
        h2->GetXaxis()->SetRangeUser(800, 1000);
      }
      else {
        h2->GetXaxis()->SetRangeUser(1400, 2300);
      }
      h2->Write();
      //outFile_FitHistoAlpha->Close();


      return true ;
    }
}


/////////////////////////////////
Double_t Calib_ZeroForceMethod(string xy,TGraphErrors *gr,float Pedestal, Double_t *mean, Double_t *sigma)
{
  // not adapted for W1

	Double_t energy[3];
	Double_t errors[3];

	if(xy=="FRONT")
		for(int i = 0 ; i < 3 ; i ++)
		{
			energy[i] = energyFront[i];
			errors[i] = errorsFront[i];
		}

	if(xy=="BACK")
		for(int i = 0 ; i < 3 ; i ++)
		{
			energy[i] = energyBack[i];
			errors[i] = errorsBack[i];
		}

	gr->SetPoint(0,Pedestal,energy[0]);

	for (Int_t p = 0; p < 3; p++) {
		gr->SetPoint(p, mean[p], energy[p]);
		gr->SetPointError(p, sigma[p], errors[p]);    
	}

	TF1 *f1 = new TF1("f1",Form("[0]*(x-%f)",Pedestal));
	gr->Fit("f1", "Q" );

	a = f1 -> GetParameter(0);

	if (xy == "X")
		calib_file << "MUST2_T" << CurrentTelescope << "_Si_X" << CurrentStrip << "_E " << b << " " << a  << endl ;

	else if (xy == "Y")
		calib_file << "MUST2_T" << CurrentTelescope << "_Si_Y" << CurrentStrip << "_E " << b << " " << a  << endl ;

	delete f1;
	return a ;
}

/////////////////////////////////
Double_t Calib_ZeroExtrapolationMethod(TH1F* hist , string xy, TGraphErrors *gr, float Pedestal, Double_t* mean, Double_t* error_mean, Double_t* sigma, Double_t &a , Double_t &b)
{  
	Double_t energy[3];
	Double_t errors[3];

	if(xy=="FRONT")
		for(int i = 0 ; i < 3 ; i ++)
		{
			energy[i] = energyFront[i];
			errors[i] = errorsFront[i];
		}

	if(xy=="BACK")
		for(int i = 0 ; i < 3 ; i ++)
		{
			energy[i] = energyBack[i];
			errors[i] = errorsBack[i];
		}

	for (Int_t p = 0; p < 3; p++) {
		gr->SetPoint(p, mean[p], energy[p]);
		gr->SetPointError(p, error_mean[p], errors[p]); // case error for mean = error fit on mean
		//gr->SetPointError(p, sigma[p], errors[p]); // case error for mean = sigma
	}

	TF1 *f1 = new TF1("f1","[1]+[0]*x");
  TFitResultPtr r = gr->Fit("f1", "S");
  Int_t fitStatus = r;
  cout << "status " << fitStatus << endl;

	a = f1 -> GetParameter(0);
	b = f1 -> GetParameter(1);
  //cout << "a = " << a << " b = " << b << endl;

	if(RefitWithSatellite)
	{
		Find_Satellites(hist);

		for (Int_t p = 0; p < 3; p++) 
		{
			gr->SetPoint(p, mean[p], energy[p]);
			gr->SetPointError(p, sigma[p], a*sigma[p]);    
      cout << "mean alpha refit satellites = " << mean[p] << " energy theo = " << energy[p] << endl;
		}

		gr->Fit("f1", "Q" );

		a = f1 -> GetParameter(0);
		b = f1 -> GetParameter(1);
	}
  //cout << "a refit = " << a << " b refit = " << b << endl;

  // save graph in a root file
  TFile *outFile_calib = new TFile("/projet/astronuc/31P3Het/nptool/Projects/nsi91/W1_Analysis/raw/graphs_calib.root", "update");
//  TFile *outFile_calib = new TFile("/projet/astronuc/31P3Het/nptool/Projects/nsi91/W1_Analysis/raw/graphs_calib.root", "recreate");
  gr->SetNameTitle(Form("grCalib_D%d_strip0%d", CurrentTelescope, CurrentStrip), Form("D%d_Strip0%d", CurrentTelescope, CurrentStrip));
  gr->Write();
  outFile_calib->Close();

	//if( (a < 9.0 &&  a > 2.0) || (a > -9.0 &&  a < -2.0) )  
	//if( (a < 0.009 &&  a > 0.006) || (a > -0.009 &&  a < -0.006) )  
		coeff_a->SetPoint(CurrentStrip,CurrentStrip,a);

	//if( (b < -54 && b > -72) || (b > 54 && b < 72) )
		coeff_b->SetPoint(CurrentStrip,CurrentStrip,b);

	// look at the dispersion around Pedestals
	Double_t dispersion = Pedestal + b/a ;
  //cout << "ADC(0) = " << -b/a << " dispersion = " << dispersion << endl;
  cout << endl;
	dispersion_file  << "W1_D" << CurrentTelescope << "_Si_F" << CurrentStrip << "_E_Zero_Dispersion " << dispersion << endl ;

	// Condition avoid Mean problem due to a few large value
	//if(dispersion< 150 && dispersion> -150 )
	//if(dispersion<60 && dispersion>-60 )
	if(dispersion<50 && dispersion>-50 )
		ZeroDispersion->SetPoint(CurrentStrip,CurrentStrip,dispersion);

	Dispersion->Fill(dispersion);

	//if(dispersion > 150 || dispersion < -150)
	//if(dispersion > 60 || dispersion < -60)
	if(dispersion > 50 || dispersion < -50)
	{
		ostringstream disp;
		disp << dispersion ;
		BadStrip[CurrentStrip] += " zero extrapolation too high:" + disp.str() +"channels; ";
	}


	if (xy == "FRONT") {
		calib_file << "W1_D" << CurrentTelescope << "_FRONT" << CurrentStrip-1 << "_E " << b << " " << a  << endl ;
    calibError_file << "W1_D" << CurrentTelescope << "_FRONT" << CurrentStrip-1 << "_E " << b << " " << f1 -> GetParError(1) <<  " " << a << " " << f1 -> GetParError(0) << endl;
  }

  else if (xy == "BACK") {
		calib_file << "W1_D" << CurrentTelescope << "_BACK" << CurrentStrip-1 << "_E " << b << " " << a  << endl ;
    calibError_file << "W1_D" << CurrentTelescope << "_BACK" << CurrentStrip-1 << "_E " << b << " " << f1 -> GetParError(1) <<  " " << a << " " << f1 -> GetParError(0) << endl;
  }

	delete f1;
	return a ;

}

/////////////////////////////////////////
void LatexSummaryHeader(TString xy)
{

	latex_file.open(folder+"/latex/"+main_name+".tex");

	///// Write File Header

	latex_file << "\\documentclass[a4paper,6pt]{article}" << endl ;
	latex_file << "\\usepackage[french]{babel}" << endl ;
	latex_file << "\\usepackage[T1]{fontenc}" << endl ;
	latex_file << "\\usepackage{graphicx}" << endl ;
	latex_file << "\\usepackage{fullpage}" << endl ;
	latex_file << "\\topmargin = 0pt" << endl ;
	latex_file << "\\headsep = 0pt" << endl ;

	// Start Document
	latex_file << "\\begin{document}" << endl ;
	latex_file << "\\title{W1 DSSD Energy Calibration Report}" << endl ;
	latex_file << "\\date{}" << endl ;
	latex_file << "\\maketitle" << endl ;

	// Write Report header
	latex_file << "\\section{Calibration Summary}" << endl ;
	latex_file << "\\begin{itemize}" << endl ;
	latex_file << "\t \\item[{\\bf Experiment:}] "<< Experiment << endl ;
	latex_file << "\t \\item[{\\bf Operator:}] "<< Operator << endl ;
	latex_file << "\t \\item[{\\bf App. Date:}] "<< Run_Period << endl ;
	latex_file << "\t \\item[{\\bf Source:}] "<< Source << endl ;
	latex_file << "\t \\item[{\\bf Dead Layer:}] "<< "Al "<< AlThickness/micrometer << "$\\mu$m + Si " << SiThickness/micrometer << "$\\mu$m" << endl ;
	latex_file << "\t \\item[{\\bf Comment:}] "<< Comment << endl ;
	latex_file << "\t \\item[] "<< endl ;
	latex_file << "\t \\item[{\\bf Calibration Method:}] "<< " " << method << " "<< endl ;
	latex_file << "\t \\item[{\\bf Telescope Treated:}] "<<  " " << Telescope_Number << endl ;
	latex_file << "\t \\item[{\\bf Strip Treated:}] "<<  " " << Strip_Start << " to "<< Strip_End << " " << endl ;
	latex_file << "\t \\item[{\\bf DSSD Side:}] "<< " " <<  xy << endl ;

	latex_file << "\\end{itemize}" << endl ;

	latex_file << "\\begin{itemize}" << endl ;
	latex_file << "\t \\item[] "<< endl ;
	latex_file << "\t \\item[] "<< endl ;
	latex_file << "\\end{itemize}" << endl ;

	latex_file << "{\\bf Source Description:} " << endl ;
	latex_file << "\\begin{center}"<<endl ;
	latex_file << "\\begin{tabular}{ | c | c | c | } "<<endl ;
	latex_file << "\\hline "<<endl ;
	latex_file << "Isotope & Original Energy (MeV) & Branching Ratio \\\\ \\hline " << endl ;

	for(int hh = 0 ; hh < Source_Number_Peak ; hh++)
	{
		latex_file << Source_isotope[hh] << " & " << Source_E[hh] << " & " << Source_branching_ratio[hh] << " \\\\ \\hline" << endl;
	}

	latex_file << "\\end{tabular} "<<endl ;
	latex_file << "\\end{center}"<<endl ;

	latex_file <<"\\pagebreak"<<endl ;
}

///
void LatexSummaryEnder()
{
	latex_file << endl <<  "\\end{document}" << endl ;
	latex_file.close();
	// generate the pdf file and clean-up
	int returnValue7 = system("pdflatex "+folder+"/latex/"+main_name+".tex");
	int returnValue8 = system("rm -f *.log");
	int returnValue9 = system("rm -f *.aux");
	int returnValue10 = system("mv " + main_name+".pdf "+folder  );
}

///
void LatexSummaryTelescope()
{
	/// Write main summary
	latex_file << "\\section{Telescope "<< CurrentTelescope << " }"<<endl ;
	/// List symptomatic strips and reason

	if(BadStrip.size()>0)
	{
		latex_file << " Bad Strip:" << endl ;
		latex_file << "\\begin{center}"<<endl ;
		latex_file << "\\begin{tabular}{ | c | c | } "<<endl ;
		latex_file << "\\hline "<<endl ;
		latex_file << " Strip Number & Problem \\\\ \\hline "<<endl ;
		map<int,string>::iterator it ;
		for(it = BadStrip.begin() ; it!=BadStrip.end() ; it++)
		{
			latex_file << it->first << " & " << it->second <<  " \\\\ \\hline "<<endl ;
		}

		latex_file << "\\end{tabular} "<<endl ;
		latex_file << "\\end{center}"<<endl ;
	}

	else
		latex_file << "Bad Strip : All Strip are ok."<<endl ;

	// Add the Graph
	TString filename = Tsummary->GetName();
	TString path = folder+"/latex/pictures/"+filename+".pdf";

	latex_file <<"\\begin{figure}[htcb!]"<<endl ;
	latex_file <<"\\begin{center}"<<endl ;
	latex_file <<"\\includegraphics[width=0.7\\textwidth]{"+path +"}"<<endl ;
	latex_file <<"\\end{center}"<<endl ;
	latex_file <<"\\end{figure}"<<endl ;

	latex_file <<"\\pagebreak"<<endl ;

	/// add summary graph and image

}


//////// Satellite finder and description of the Peak+Sattelite look-a-like function
void Find_Satellites(TH1F *h)
{

	if(mean[0]==0 && mean[1]==0 && mean[2]==0) { cout << "pas de pics ---> pas de satellites!" << endl;}

	else {

		Float_t linf1 =0 , lsup1 =0, linf2 =0 , lsup2 =0 , linf3 =0 , lsup3=0;

		if(a>0) { // ie Y case
			//linf1 = mean[0]-15; lsup1 = mean[0]+10; // adjust ?
/*			linf1 = mean[0]-15; lsup1 = mean[0]+10; // adjust ?
			linf2 = mean[1]-15; lsup2 = mean[1]+10;
			linf3 = mean[2]-15; lsup3 = mean[2]+10;*/
			linf1 = mean[0]-25; lsup1 = mean[0]+10; // adjust ?
			linf2 = mean[1]-25; lsup2 = mean[1]+10;
			linf3 = mean[2]-25; lsup3 = mean[2]+10;
		}

		else { // ie X case 
			lsup1 = mean[0]+15; linf1 = mean[0]-10;
			lsup2 = mean[1]+15; linf2 = mean[1]-10;
			lsup3 = mean[2]+15; linf3 = mean[2]-10;
		}

		Double_t keVtoMeV = 1./1000. ;

		TF1 *Pu = new TF1("fit_sat_Pu", source_Pu, linf1, lsup1, 6);
		Pu->SetParameters(150,mean[0],mean[0]-12.4*keVtoMeV/a,mean[0]-51.6*keVtoMeV/a,sigma[0]);
		Pu->SetParLimits(2,mean[0]-12.4*keVtoMeV/a-10,mean[0]-12.6*keVtoMeV/a+10);
		Pu->SetParLimits(3,mean[0]-51.6*keVtoMeV/a-10,mean[0]-51.6*keVtoMeV/a+10);
		Pu->SetParNames("Constant","Mean_value1","Mean_value2","Mean_value3","SigmaPu");
		h->Fit("fit_sat_Pu", "RQ");

		TF1 *Am = new TF1("fit_sat_Am", source_Am, linf2, lsup2, 6);
		Am->SetParameters(150,mean[1],mean[1]-43.2*keVtoMeV/a,mean[1]-98.4*keVtoMeV/a,sigma[1]);
		Am->SetParLimits(2,mean[1]-43.2*keVtoMeV/a-10,mean[1]-43.2*keVtoMeV/a+10);
		Am->SetParLimits(3,mean[1]-98.4*keVtoMeV/a-10,mean[1]-98.4*keVtoMeV/a+10);
		Am->SetParNames("Constant","Mean_value1","Mean_value2","Mean_value3","SigmaAm");
		h->Fit("fit_sat_Am", "RQ+");

		TF1 *Cm = new TF1("fit_sat_Cm", source_Cm, linf3, lsup3, 6);
		Cm->SetParameters(150,mean[2],mean[2]-43.1*keVtoMeV/a,sigma[2]);
		Cm->SetParLimits(2,mean[2]-43.1*keVtoMeV/a-10,mean[0]-43.1*keVtoMeV/a-10);
		Cm->SetParNames("Constant","Mean_value1","Mean_value2","SigmaCm");
		h->Fit("fit_sat_Cm", "RQ+");

		mean[0]=Pu->GetParameter(1);  // Position of the 1st principal peak
		sigma[0]=Pu->GetParameter(4); // Sigma of the 1st principal peak
		sigma_fit->Fill(sigma[0]) ;
		error_par[0]= Pu->GetParError(1);
		mean[1]=Am->GetParameter(1);
		sigma[1]=Am->GetParameter(4);
		sigma_fit->Fill(sigma[1]) ;
		error_par[1]= Am->GetParError(1);
		mean[2]=Cm->GetParameter(1);
		sigma[2]=Cm->GetParameter(3);
		sigma_fit->Fill(sigma[2]) ;
		error_par[2]= Cm->GetParError(1);

    // save graph and fit in a root file
    TFile *outFile_FitSatellite_HistoAlpha = new TFile("/projet/astronuc/31P3Het/nptool/Projects/nsi91/W1_Analysis/raw/graphs_FitSatelliteAlpha.root", "update");
    h->SetNameTitle(Form("grFitSatelliteAlpha_D%d_Strip0%d", CurrentTelescope, CurrentStrip), Form("D%d_Strip0%d", CurrentTelescope, CurrentStrip));
    h->Write();
    outFile_FitSatellite_HistoAlpha->Close();

	}

}

///////////////////////////////////////////////
Double_t source_Pu(Double_t *x, Double_t *par)
{
	// [0] : constant
	// [1] : position peak1
	// [2] : position peak2
	// [3] : position peak3
	// [4] : sigma

	Double_t arg1 = 0;
	Double_t arg2 = 0;
	Double_t arg3 = 0;

	if(par[4]!=0) { 
		arg1 = (x[0]-par[1])/par[4];
		arg2 = (x[0]-par[2])/par[4];
		arg3 = (x[0]-par[3])/par[4];
	}

	else cout << " Attention, sigma est nul !" << endl;

	Double_t gaus1 =           par[0]*exp(-0.5*arg1*arg1);
	Double_t gaus2 = 15.1/73.8*par[0]*exp(-0.5*arg2*arg2);
	Double_t gaus3 = 11.5/73.8*par[0]*exp(-0.5*arg3*arg3);
	Double_t fitval = gaus1+gaus2+gaus3;

	return fitval;
}

///////////////////////////////////////////////
Double_t source_Am(Double_t *x, Double_t *par)
{
	// [0] : constant
	// [1] : position peak1
	// [2] : position peak2
	// [3] : position peak3
	// [4] : sigma

	Double_t arg1 = 0;
	Double_t arg2 = 0;
	Double_t arg3 = 0;

	if(par[4]!=0) { 
		arg1 = (x[0]-par[1])/par[4];
		arg2 = (x[0]-par[2])/par[4];
		arg3 = (x[0]-par[3])/par[4];
	}

	else cout << " Attention, sigma est nul !" << endl;

	Double_t gaus1 =           par[0]*exp(-0.5*arg1*arg1);
	Double_t gaus2 = 13.0/84.5*par[0]*exp(-0.5*arg2*arg2);
	Double_t gaus3 = 1.6/84.5 *par[0]*exp(-0.5*arg3*arg3);
	Double_t fitval= gaus1+gaus2+gaus3;

	return fitval;
}

///////////////////////////////////////////////
Double_t source_Cm(Double_t *x, Double_t *par)
{
	// [0] : constante
	// [1] : position peak1
	// [2] : position peak2
	// [3] : sigma

	Double_t arg1 = 0;
	Double_t arg2 = 0;

	if(par[3]!=0) { 
		arg1 = (x[0]-par[1])/par[3];
		arg2 = (x[0]-par[2])/par[3];
	}

	else cout << " Attention, sigma est nul !" << endl;

	Double_t gaus1 =           par[0]*exp(-0.5*arg1*arg1);
	Double_t gaus2 = 23.6/76.4*par[0]*exp(-0.5*arg2*arg2);
	Double_t fitval= gaus1+gaus2; 

	return fitval;
}  
