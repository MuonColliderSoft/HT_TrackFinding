/*
	
			
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////			
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

This macro is supposed to produce all plots needed for the paper in their final format.
Input is taken from a number of .root files containing histograms and graphs.
Output is produced in .png form and stored in the "plots" folder. An individual .png file 
is produced for each plot. An Index.txt file is also produced with the list of the files 
with their contents.

// to run under root:

//  root> .x makeAllPlots.C
				
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////			
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


*/

	// Histogram files //////////////////////////////////////////////////////
		
		string histFileNameFullEff = "AAAEventGeneration10KFullEffEta2.5.root";
		string histFileName1pcIneff = "AAAEventGeneration10K1pcIneffEta2.5.root";		
		string histFileNameFullEffEta2 = "AAAEventGeneration10KFullEffEta2.root";
		string histFileName1pcIneffEta2 = "AAAEventGeneration10K1pcIneffEta2.root";
		string histFileName60ps = "AAAAParticleIDTests60ps.root";
		string histFileName10ps = "AAAAParticleIDTests10ps.root";
		string histFileName1ps = "AAAAParticleIDTests1ps.root";
		string histFileNameChi2 = "AAA10KEventGeneration.root";
		string histFileNameHiStatChi2 = "AAAEventGeneration_15_360_6-Pt_3p0-HitIneff_0p01-xphiErr_10um-2M.root";
		
	/////////////////////////////////////////////////////////////////////////		

	// canvasses are automatically numbered starting from 1
		
	static int indCanvas = 0;
	
	ofstream indexFile; // this is the index.txt file
	
	
	
//////////////////////////////////////////////////////////////////////////////////////	
//////////////////////////////////////////////////////////////////////////////////////
// A call to getCanvas will create a new canvas and provide a pointer to it (canvasStar)
	
	void getCanvas (TCanvas * &canvasStar){	
	
		// the canvas will be shown on the screen, each one at a different position

		// first canvas position and size

		static int canv_x = 50;
		static int canv_y = 50;
		static int canv_width = 500;
		static int canv_height = 400;
		
		int verticalShift = 40; // vertical increment of canvas position
		
		canv_y += verticalShift;
				
		++indCanvas;
		string canvasName = "c"+ std::to_string(indCanvas); //synthesize canvas name
				
		// create new canvas
		canvasStar = new TCanvas(canvasName.c_str(),canvasName.c_str(),canv_x,canv_y,canv_width,canv_height);
		
		return;
	
	}
	
	
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// A call to getGraph will search for a graph by name (plotname) in a .root file 
// (histFile) and return a pointer to it (graphStar).	 
	
	void getGraph (TFile *histFile, string plotName, TGraph * &graphStar){	

		graphStar = (TGraph*)histFile->Get(plotName.c_str());
		if(graphStar == NULL) cout << "ERROR *** graph " << plotName << " not found" << endl;
		
		return;
		
	}
		
	
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// A call to getHist will search for a histogram by name (plotname) in a .root file 
// (histFile) and return a pointer to it (histStar).
	
	void getHist(TFile *histFile, string plotName, TH1D * &histStar){
	
		histStar =(TH1D*)histFile->Get(plotName.c_str());
		if(histStar == NULL) cout << "ERROR *** hist " << plotName << " not found" << endl;
		
		return;
	
	}
	
	
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// setStyleEffGraph does all the necessary style adjustments for efficiency plots	

	void setStyleEffGraph (TGraph * graphStar, string plotName, string plotTitle, TCanvas * canvasStar) {
	
		canvasStar->SetGrid();
		
		graphStar->SetTitle("");				
		graphStar->SetMarkerStyle(8);						
		graphStar->SetMarkerSize(1.5);
		graphStar->SetMaximum(1.);
		graphStar->SetMinimum(0.95);
		
		return;		
	}
	
		
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// setStyleEffGraphRed does all the necessary style adjustments for efficiency plots	
// to be superimposed with hollow circle markers (4)
	
	void setStyleEffGraphSup (TGraph * graphStar, string plotName, string plotTitle, TCanvas * canvasStar) {
		
		canvasStar->SetGrid();
		
		graphStar->SetTitle("");				
		graphStar->SetMarkerStyle(4);						
		graphStar->SetMarkerSize(1.5);
		graphStar->SetMaximum(1.);
		
		return;		
	}
	
	
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// setStyleResHist does all the necessary style adjustments for resolution histograms	
		
	void setStyleResHist (TH1D * histStar, string plotName, string plotTitle, TCanvas * canvasStar) {

		canvasStar->SetGrid();
		
		histStar->SetTitle("");	
		histStar->SetMarkerStyle(8);						
		histStar->SetMarkerSize(1.5);
		histStar->SetLineColor(1);
		histStar->SetMinimum(0.);
		
		return;		
	}
	
	
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// setStyleDistHist does all the necessary style adjustments for distribution histograms	
		
	void setStyleDistHist (TH1D * histStar, string plotName, string plotTitle, TCanvas * canvasStar) {
		
		histStar->SetTitle("");	
		histStar->SetLineWidth(2);
		histStar->GetXaxis()->CenterTitle(true);
		histStar->GetYaxis()->CenterTitle(true);
		histStar->SetStats(0);	
		histStar->SetMarkerStyle(8);						
		histStar->SetMarkerSize(1.5);
		histStar->SetLineColor(1);
		histStar->SetMinimum(0.);
		
		return;		
	}

		
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// PrintGraph creates a .png file with the plot pointed by graphStar. 
// The .png file is stored in the "plots" folder. A line is also added to the index.txt file.
	
	void PrintGraph (TGraph * graphStar, string plotName, string plotTitle, TCanvas * canvasStar) {
			
		string printName = "plots/"+plotName+".png";		
		canvasStar->SaveAs(printName.c_str());
		cout << indCanvas << ". " << plotTitle << endl;
		indexFile << indCanvas << ". " << plotName+".png : " << plotTitle << endl;
		
		return;	
	}
	
			
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// PrintHist creates a .pdf file with the histogram pointed by histStar. 
// The .pdf file is stored in the "plots" folder. A line is also added to the index.txt file.
		
	void PrintHist (TH1D * histStar, string plotName, string plotTitle, TCanvas * canvasStar) {
			
		string printName = "plots/"+plotName+".png";		
		canvasStar->SaveAs(printName.c_str());
		cout << indCanvas << ". " << plotTitle << endl;
		indexFile << indCanvas << ". " << plotName+".png : " << plotTitle << endl;
		
		return;	
	}

			
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////			
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////			
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////			
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////			
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////			
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
	
	
	
	void makeAllPlots(){
	
		gErrorIgnoreLevel = kWarning;
	
		
		TCanvas * canvasStar;
		TGraph * graphStar;
		TGraph * graphStar1pcIneff;
		TH1D * histStar;
		TH1D * histStar2;
		string plotName;
		string plotTitle;
		
		
  		indexFile.open ("plots/Index.txt");
  			indexFile << "------ Plot file index ------\n\n";
  	
 
		TFile *histFileFullEff = new TFile(histFileNameFullEff.c_str());
		if (!histFileFullEff) {cout << " histogram file " << histFileNameFullEff << " not found" << endl;
						return;
						}
			
		TFile *histFile1pcIneff = new TFile(histFileName1pcIneff.c_str());
		if (!histFile1pcIneff) {cout << " histogram file " << histFileName1pcIneff << " not found" << endl;
						return;
						}
	
		TFile *histFileFullEffEta2 = new TFile(histFileNameFullEffEta2.c_str());
		if (!histFileFullEffEta2) {cout << " histogram file " << histFileNameFullEffEta2 << " not found" << endl;
						return;
						}
			
		TFile *histFile1pcIneffEta2 = new TFile(histFileName1pcIneffEta2.c_str());
		if (!histFile1pcIneffEta2) {cout << " histogram file " << histFileName1pcIneffEta2 << " not found" << endl;
						return;
						}
						
		
//////////////////////////////////////////////////////////////////////////////////////
//
//	NUMBER OF HITS VS. ETA
//
//////////////////////////////////////////////////////////////////////////////////////
	
		plotName = "HTrackNhitsVsEta";
		plotTitle = "Number of hits vs eta";
					
		
			getCanvas(canvasStar);
			getHist(histFileFullEff, plotName, histStar);
			setStyleResHist (histStar, plotName, plotTitle, canvasStar);			
			histStar->GetYaxis()->SetRangeUser(3.,9.);
			histStar->GetXaxis()->SetTitle("eta");
			histStar->GetYaxis()->SetTitle("Nhits");	
			histStar->GetXaxis()->CenterTitle(true);
			histStar->GetYaxis()->CenterTitle(true);			
			histStar->SetStats(0);
			histStar->Draw(); 
								
			PrintHist (histStar, plotName, plotTitle, canvasStar);
			
	
//////////////////////////////////////////////////////////////////////////////////////
//
//	EFFICIENCY VS. ETA
//
//////////////////////////////////////////////////////////////////////////////////////

	// This is actually the superposition of two plots: one with single hit efficiency
	// of 100% (from histFileFullEff) and one with single hit efficiency of 99% (from
	// histFile1pcIneff)
		
		plotName = "effVsEta";
		plotTitle = "Efficiency vs eta";
					
			getCanvas(canvasStar);
			getGraph(histFileFullEff, plotName, graphStar);	
			setStyleEffGraph(graphStar, plotName, plotTitle, canvasStar);
			graphStar->GetXaxis()->SetTitle("eta");
			graphStar->GetYaxis()->SetTitle("efficiency");	
			graphStar->GetXaxis()->CenterTitle(true);
			graphStar->GetYaxis()->CenterTitle(true);
			graphStar->Draw("AP"); 
			
			getGraph(histFile1pcIneff, plotName, graphStar1pcIneff);	
			setStyleEffGraphSup (graphStar1pcIneff, plotName, plotTitle, canvasStar);				
			graphStar1pcIneff->Draw("PSAME");
			
			auto legend = new TLegend(0.2,0.2);
   			legend->AddEntry(graphStar,"Single hit full eff","lep");
   			legend->AddEntry(graphStar1pcIneff,"Single hit 1% ineff","lep");
   			legend->Draw();	
   			
   			PrintGraph (graphStar1pcIneff, plotName, plotTitle, canvasStar);
		
		
//////////////////////////////////////////////////////////////////////////////////////
//
//	EFFICIENCY VS. 1/PT (eta range = [-2, +2])
//
//////////////////////////////////////////////////////////////////////////////////////

	// This is actually the superposition of two plots: one with single hit efficiency
	// of 100% (from histFileFullEff) and one with single hit efficiency of 99% (from
	// histFile1pcIneff)
		
		plotName = "effVsInvPt";
		plotTitle = "Efficiency vs 1/PT";
			
			getCanvas(canvasStar);
			getGraph(histFileFullEffEta2, plotName, graphStar);				
			setStyleEffGraph (graphStar, plotName, plotTitle, canvasStar);	
			graphStar->GetXaxis()->SetTitle("1/PT (GeV/c-1)");
			graphStar->GetYaxis()->SetTitle("efficiency");	
			graphStar->GetXaxis()->CenterTitle(true);
			graphStar->GetYaxis()->CenterTitle(true);
			graphStar->Draw("AP"); 
			
			getGraph(histFile1pcIneffEta2, plotName, graphStar1pcIneff);	
			setStyleEffGraphSup (graphStar1pcIneff, plotName, plotTitle, canvasStar);						
			graphStar1pcIneff->Draw("PSAME");
					
			legend = new TLegend(0.2,0.2);
   			legend->AddEntry(graphStar,"Single hit full eff","lep");
   			legend->AddEntry(graphStar1pcIneff,"Single hit 1% ineff","lep");
   			legend->Draw();
   			
   			PrintGraph (graphStar1pcIneff, plotName, plotTitle, canvasStar);
  
//////////////////////////////////////////////////////////////////////////////////////
//
//	PT RESOLUTION VS. ETA
//
//////////////////////////////////////////////////////////////////////////////////////
								
		plotName = "HSigmaInvPtVsEta";
		plotTitle = "sigmaPT/PT^2 Vs. Eta";
		
			getCanvas(canvasStar);
			getHist(histFileFullEff, plotName, histStar);
			setStyleResHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetTitle("eta");
			histStar->GetYaxis()->SetTitle("sigmaPT/PT^2 (GeV/c-1)");	
			histStar->GetXaxis()->CenterTitle(true);
			histStar->GetYaxis()->CenterTitle(true);
			histStar->Draw(); 
								
			PrintHist (histStar, plotName, plotTitle, canvasStar);
			
//////////////////////////////////////////////////////////////////////////////////////
//
//	PT RESOLUTION VS. 1/PT (eta range = [-2, +2])
//
//////////////////////////////////////////////////////////////////////////////////////
											
		plotName = "HSigmaInvPtVsInvPt";
		plotTitle = "sigmaPT/PT^2 Vs. 1/PT";
		
			getCanvas(canvasStar);
			getHist(histFileFullEffEta2, plotName, histStar);
			setStyleResHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetTitle("1/PT (GeV/c-1)");
			histStar->GetYaxis()->SetTitle("sigmaPT/PT^2 (GeV/c-1)");	
			histStar->GetXaxis()->CenterTitle(true);
			histStar->GetYaxis()->CenterTitle(true);
			histStar->Draw();
			 				
			PrintHist (histStar, plotName, plotTitle, canvasStar);
			
			
						
		
	
//////////////////////////////////////////////////////////////////////////////////////
//
//	EFFICIENCY VS. ABS(ETA)
//
//////////////////////////////////////////////////////////////////////////////////////

	// This is actually the superposition of two plots: one with single hit efficiency
	// of 100% (from histFileFullEff) and one with single hit efficiency of 99% (from
	// histFile1pcIneff)
		
		plotName = "effVsAbsEta";
		plotTitle = "Efficiency vs abs(eta)";
					
			getCanvas(canvasStar);
			getGraph(histFileFullEff, plotName, graphStar);	
			setStyleEffGraph(graphStar, plotName, plotTitle, canvasStar);
			graphStar->GetXaxis()->SetTitle("eta");
			graphStar->GetYaxis()->SetTitle("efficiency");	
			graphStar->GetXaxis()->CenterTitle(true);
			graphStar->GetYaxis()->CenterTitle(true);
			graphStar->Draw("AP"); 
			
			getGraph(histFile1pcIneff, plotName, graphStar1pcIneff);	
			setStyleEffGraphSup (graphStar1pcIneff, plotName, plotTitle, canvasStar);				
			graphStar1pcIneff->Draw("PSAME");
			
			auto legend2 = new TLegend(0.2,0.2);
   			legend2->AddEntry(graphStar,"Single hit full eff","lep");
   			legend2->AddEntry(graphStar1pcIneff,"Single hit 1% ineff","lep");
   			legend2->Draw();	
   			
   			PrintGraph (graphStar1pcIneff, plotName, plotTitle, canvasStar);
		
		
//////////////////////////////////////////////////////////////////////////////////////
//
//	EFFICIENCY VS. ABS(1/PT) (eta range = [-2, +2])
//
//////////////////////////////////////////////////////////////////////////////////////

	// This is actually the superposition of two plots: one with single hit efficiency
	// of 100% (from histFileFullEff) and one with single hit efficiency of 99% (from
	// histFile1pcIneff)
		
		plotName = "effVsAbsInvPt";
		plotTitle = "Efficiency vs abs(1/PT)";
			
			getCanvas(canvasStar);
			getGraph(histFileFullEffEta2, plotName, graphStar);				
			setStyleEffGraph (graphStar, plotName, plotTitle, canvasStar);	
			graphStar->GetXaxis()->SetTitle("1/PT (GeV/c-1)");
			graphStar->GetYaxis()->SetTitle("efficiency");	
			graphStar->GetXaxis()->CenterTitle(true);
			graphStar->GetYaxis()->CenterTitle(true);
			graphStar->Draw("AP"); 
			
			getGraph(histFile1pcIneffEta2, plotName, graphStar1pcIneff);	
			setStyleEffGraphSup (graphStar1pcIneff, plotName, plotTitle, canvasStar);						
			graphStar1pcIneff->Draw("PSAME");
					
			legend = new TLegend(0.2,0.2);
   			legend->AddEntry(graphStar,"Single hit full eff","lep");
   			legend->AddEntry(graphStar1pcIneff,"Single hit 1% ineff","lep");
   			legend->Draw();
   			
   			PrintGraph (graphStar1pcIneff, plotName, plotTitle, canvasStar);
  
//////////////////////////////////////////////////////////////////////////////////////
//
//	PT RESOLUTION VS. ABS(ETA)
//
//////////////////////////////////////////////////////////////////////////////////////
								
		plotName = "HSigmaInvPtVsAbsEta";
		plotTitle = "sigmaPT/PT^2 Vs. abs(Eta)";
		
			getCanvas(canvasStar);
			getHist(histFileFullEff, plotName, histStar);
			setStyleResHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetTitle("eta");
			histStar->GetYaxis()->SetTitle("sigmaPT/PT^2 (GeV/c-1)");	
			histStar->GetXaxis()->CenterTitle(true);
			histStar->GetYaxis()->CenterTitle(true);
			histStar->Draw(); 
								
			PrintHist (histStar, plotName, plotTitle, canvasStar);
			
//////////////////////////////////////////////////////////////////////////////////////
//
//	PT RESOLUTION VS. ABS(1/PT) (eta range = [-2, +2])
//
//////////////////////////////////////////////////////////////////////////////////////
											
		plotName = "HSigmaInvPtVsAbsInvPt";
		plotTitle = "sigmaPT/PT^2 Vs. abs(1/PT)";
		
			getCanvas(canvasStar);
			getHist(histFileFullEffEta2, plotName, histStar);
			setStyleResHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetTitle("1/PT (GeV/c-1)");
			histStar->GetYaxis()->SetTitle("sigmaPT/PT^2 (GeV/c-1)");	
			histStar->GetXaxis()->CenterTitle(true);
			histStar->GetYaxis()->CenterTitle(true);
			histStar->Draw();
			 				
			PrintHist (histStar, plotName, plotTitle, canvasStar);
		

		
			
//////////////////////////////////////////////////////////////////////////////////////
//
//	TRACK PARAMETER DISTRIBUTIONS
//
//////////////////////////////////////////////////////////////////////////////////////
		
		plotName = "HTrackZ0";
		plotTitle = "Track origin z position";
				
			getCanvas(canvasStar);
			getHist(histFileFullEff, plotName, histStar);
			setStyleDistHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetRangeUser(-20,+20);
			histStar->GetXaxis()->SetTitle("Z (mm)");
			histStar->SetStats(1);
			gStyle->SetOptStat("r");
			histStar->Draw();
			 				
			PrintHist (histStar, plotName, plotTitle, canvasStar);
			
		
		plotName = "HTrackT0";
		plotTitle = "Track origin time";
					
			getCanvas(canvasStar);
			getHist(histFileFullEff, plotName, histStar);
			setStyleDistHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetRangeUser(-20,+20);
			histStar->GetXaxis()->SetTitle("T0 (mm)");
			histStar->SetStats(1);
			gStyle->SetOptStat("r");
			histStar->Draw();
			 				
			PrintHist (histStar, plotName, plotTitle, canvasStar);
		
				
		plotName = "HTrackEta";
		plotTitle = "Track eta at origin";
					
			getCanvas(canvasStar);
			getHist(histFileFullEff, plotName, histStar);
			setStyleDistHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetRangeUser(-20,+20);
			histStar->GetXaxis()->SetTitle("eta");
			histStar->Draw();
			 				
			PrintHist (histStar, plotName, plotTitle, canvasStar);
		
						
		plotName = "HTrackPhi";
		plotTitle = "Track phi at origin";
					
			getCanvas(canvasStar);
			getHist(histFileFullEff, plotName, histStar);
			setStyleDistHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetRangeUser(-20,+20);
			histStar->GetXaxis()->SetTitle("phi");
			histStar->Draw();
			 				
			PrintHist (histStar, plotName, plotTitle, canvasStar);
		
			
		plotName = "HTrackInvPt";
		plotTitle = "Track 1/PT";
					
			getCanvas(canvasStar);
			getHist(histFileFullEff, plotName, histStar);
			setStyleDistHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetRangeUser(-20,+20);
			histStar->GetXaxis()->SetTitle("1/PT (GeV/c-1)");
			histStar->Draw();
			 				
			PrintHist (histStar, plotName, plotTitle, canvasStar);

			
//////////////////////////////////////////////////////////////////////////////////////
//
//	OPTIMIZATION OF EXECUTION SPEED
//
//////////////////////////////////////////////////////////////////////////////////////
			
			plotName = "DeltaEta";
			plotTitle = "Delta eta";
				
			getCanvas(canvasStar);
			getHist(histFileFullEff, plotName, histStar);
			setStyleDistHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetTitle("Delta eta");
			histStar->Draw();
			 				
			PrintHist (histStar, plotName, plotTitle, canvasStar);
				
			
			plotName = "DeltaPhi";
			plotTitle = "Delta phi uncorrected";
				
			getCanvas(canvasStar);
			getHist(histFileFullEff, plotName, histStar);
			setStyleDistHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetTitle("Delta phi");
			histStar->Draw();
			 				
			PrintHist (histStar, plotName, plotTitle, canvasStar);
				
				
			plotName = "DeltaPhi2";
			plotTitle = "Delta phi corrected for PT";
				
			getCanvas(canvasStar);
			getHist(histFileFullEff, plotName, histStar);
			setStyleDistHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetTitle("Delta phi");
			histStar->Draw();
			 				
			PrintHist (histStar, plotName, plotTitle, canvasStar);
				
								
							
//////////////////////////////////////////////////////////////////////////////////////
//
//	DETECTOR LAYOUT
//
//////////////////////////////////////////////////////////////////////////////////////
		
						
			plotName = "HitTRZ";
			plotTitle = "R-Z Detector Layout";
		
			getCanvas(canvasStar);
			getHist(histFileFullEff, plotName, histStar);
			histStar->SetTitle("");	
			histStar->SetMarkerStyle(8);						
			histStar->SetMarkerSize(1.5);
			//histStar->SetLineColor(1);
			//histStar->SetMinimum(0.);
			histStar->SetMarkerStyle(8);						
			histStar->SetMarkerSize(0.05);
			histStar->GetXaxis()->SetTitle("Z (mm)");
			histStar->GetYaxis()->SetTitle("R (mm)");	
			histStar->GetXaxis()->CenterTitle(true);
			histStar->GetYaxis()->CenterTitle(true);
			histStar->Draw();
			 				
			PrintHist (histStar, plotName, plotTitle, canvasStar);
			
	
								
//////////////////////////////////////////////////////////////////////////////////////
//
//	TIME VS BIB FRACTION
//
//////////////////////////////////////////////////////////////////////////////////////

							
		plotName = "ExecTime";
		plotTitle = "Execution time vs. BIB density";
  		
  		getCanvas(canvasStar);
  		canvasStar->SetGrid();	
		
		Int_t n = 10;
		Double_t bib[n], t_rec[n], t_fit[n];
		bib[0] = 0.5; t_rec[0] = 5.13; t_fit[0] = 5.12;
		bib[1] = 1.0; t_rec[1] = 10.23; t_fit[1] = 10.28;
		bib[2] = 1.5; t_rec[2] = 15.44; t_fit[2] = 15.93;
		bib[3] = 2.0; t_rec[3] = 20.51; t_fit[3] = 23.19;
		bib[4] = 2.5; t_rec[4] = 25.89; t_fit[4] = 33.80;
		bib[5] = 3.0; t_rec[5] = 31.17; t_fit[5] = 51.51;
		bib[6] = 3.5; t_rec[6] = 36.31; t_fit[6] = 80.93;
		bib[7] = 4.0; t_rec[7] = 41.84; t_fit[7] = 127.30;
		bib[8] = 4.5; t_rec[8] = 47.45; t_fit[8] = 205.11;
		bib[9] = 5.0; t_rec[9] = 52.34; t_fit[9] = 325.21;				
		
		TGraph * g_trec = new TGraph (n, bib, t_rec);
		TGraph * g_tfit = new TGraph (n, bib, t_fit);
		g_tfit->SetTitle("");
		g_tfit->SetMinimum(0.);	
		g_tfit->SetMarkerStyle(4);		
		g_tfit->GetXaxis()->SetTitle("BIB density multiplier");
		g_tfit->GetYaxis()->SetTitle("Execution time");	
		g_tfit->GetXaxis()->CenterTitle(true);
		g_tfit->GetYaxis()->CenterTitle(true);		
		g_tfit->Draw("AP"); 
		
		g_trec->SetMarkerStyle(8);
		//g_tfit->SetMarkerColor(4);
		g_trec->Draw("PSAME"); 
		
			
		legend = new TLegend(0.2,0.2);
   		legend->AddEntry(g_tfit,"With fit","p");
   		legend->AddEntry(g_trec,"Without fit","p");
   		legend->Draw();	
   		
			
		PrintGraph (g_trec, plotName, plotTitle, canvasStar);
		
				
				
//////////////////////////////////////////////////////////////////////////////////////
//
//	MASS DISTRIBUTIONS
//
//////////////////////////////////////////////////////////////////////////////////////



		TFile *histFile60ps = new TFile(histFileName60ps.c_str());
		if (!histFile60ps) {cout << " histogram file " << histFileName60ps << " not found" << endl;
						return;
						}
			
		TFile *histFile10ps = new TFile(histFileName10ps.c_str());
		if (!histFile10ps) {cout << " histogram file " << histFileName10ps << " not found" << endl;
						return;
						}
									
		TFile *histFile1ps = new TFile(histFileName1ps.c_str());
		if (!histFile10ps) {cout << " histogram file " << histFileName1ps << " not found" << endl;
						return;
						}
						
		
		plotName = "HTrackMass";
		plotTitle = "Mass Distribution 60ps";
				
			getCanvas(canvasStar);
			getHist(histFile60ps, plotName, histStar);
			setStyleDistHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetRangeUser(0.05,0.8);
			histStar->GetXaxis()->SetTitle("mass (GeV)");
			histStar->SetStats(0);
			//gStyle->SetOptStat("r");
			histStar->Draw();
			 				
			PrintHist (histStar, plotName + "60ps", plotTitle, canvasStar);
			
		
		plotName = "HTrackMass";
		plotTitle = "Mass Distribution 10ps";
				
			getCanvas(canvasStar);
			getHist(histFile10ps, plotName, histStar);
			setStyleDistHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetRangeUser(0.05,0.8);
			histStar->GetXaxis()->SetTitle("mass (GeV)");
			histStar->SetStats(0);
			//gStyle->SetOptStat("r");
			histStar->Draw();
			 				
			PrintHist (histStar, plotName + "10ps", plotTitle, canvasStar);
			
			
		plotName = "HTrackMass";
		plotTitle = "Mass Distribution 1ps";
				
			getCanvas(canvasStar);
			getHist(histFile1ps, plotName, histStar);
			setStyleDistHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetRangeUser(0.05,0.8);
			histStar->GetXaxis()->SetTitle("mass (GeV)");
			histStar->SetStats(0);
			//gStyle->SetOptStat("r");
			histStar->Draw();
			 				
			PrintHist (histStar, plotName + "1ps", plotTitle, canvasStar);
			
					
			
//////////////////////////////////////////////////////////////////////////////////////
//
//	CHI SQUARE PARAMETER DISTRIBUTIONS
//
//////////////////////////////////////////////////////////////////////////////////////

							
		TFile *histFileChi2 = new TFile(histFileNameChi2.c_str());
		if (!histFileChi2) {cout << " histogram file " << histFileNameChi2 << " not found" << endl;
						return;
						}
							
		TFile *histFileHiStatChi2 = new TFile(histFileNameHiStatChi2.c_str());
		if (!histFileHiStatChi2) {cout << " histogram file " << histFileNameHiStatChi2 << " not found" << endl;
						return;
						}
		
		plotName = "HFitChi2";
		plotTitle = "Fit Chi Squared";
				
			getCanvas(canvasStar);
			getHist(histFileChi2, plotName, histStar);
			setStyleDistHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetRangeUser(0.,100.);
			histStar->GetXaxis()->SetTitle("Chi2");
			histStar->Draw();
			getHist(histFileHiStatChi2, plotName, histStar2);
			histStar2->Draw("SAME");
			 				
			PrintHist (histStar, plotName, plotTitle, canvasStar);
					
			
			
//////////////////////////////////////////////////////////////////////////////////////

		
		
			return;
		

//////////////////////////////////////////////////////////////////////////////////////
		}
	
