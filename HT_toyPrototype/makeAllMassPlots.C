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

//  root> .x makeAllMassPlots.C
				
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////			
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


*/

	// Histogram files //////////////////////////////////////////////////////
		
		string histFileName60ps = "AAAEventGeneration_15_360_12-massFit_60ps-40k.root";
		string histFileName20ps = "AAAEventGeneration_15_360_12-massFit_20ps-40k.root";
		string histFileName10ps = "AAAEventGeneration_15_360_12-massFit_10ps-40k.root";
		string histFileName1ps = "AAAEventGeneration_15_360_12-massFit_01ps-40k.root";
		
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

		static int canv_x = 20;
		static int canv_y = 20;
		static int canv_width = 1000;
		static int canv_height = 400;
		
		int verticalShift = 0; // vertical increment of canvas position
		int horizontalShift = 20; // horizontal increment of canvas position
		
		canv_x += horizontalShift;
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
		histStar->SetStats(0);
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
	
	
	
	void makeAllMassPlots(){
	
		gErrorIgnoreLevel = kWarning;
	
		
		TCanvas * canvasStar;
		TCanvas * canvasStar2;
		TGraph * graphStar;
		TGraph * graphStar1pcIneff;
		TH1D * histStar;
		TH1D * histStar2;
		string plotName;
		string plotTitle;
		
		
  		indexFile.open ("plots/IndexMass.txt");
  			indexFile << "------ Mass plots file index ------\n\n";
				
//////////////////////////////////////////////////////////////////////////////////////
//
//	MASS DISTRIBUTIONS
//
//////////////////////////////////////////////////////////////////////////////////////



		TFile *histFile60ps = new TFile(histFileName60ps.c_str());
		if (!histFile60ps) {cout << " histogram file " << histFileName60ps << " not found" << endl;
						return;
						}
			
		TFile *histFile20ps = new TFile(histFileName20ps.c_str());
		if (!histFile20ps) {cout << " histogram file " << histFileName20ps << " not found" << endl;
						return;
						}
			
		TFile *histFile10ps = new TFile(histFileName10ps.c_str());
		if (!histFile10ps) {cout << " histogram file " << histFileName10ps << " not found" << endl;
						return;
						}
									
		TFile *histFile1ps = new TFile(histFileName1ps.c_str());
		if (!histFile1ps) {cout << " histogram file " << histFileName1ps << " not found" << endl;
						return;
						}
						
		
		plotName = "HTrackMass";
		plotTitle = "Mass Distribution 60ps";
				
			getCanvas(canvasStar);
			getHist(histFile60ps, plotName, histStar);
			setStyleDistHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetRangeUser(0.01,0.8);
			histStar->GetXaxis()->SetTitle("mass (GeV)");
			histStar->SetStats(0);
			//gStyle->SetOptStat("r");
			histStar->Draw();
			
			auto legend = new TLegend(0.1,0.8,0.3,0.9);			
   			legend->AddEntry((TObject*)0, "Time resolution 60ps", "");
   			legend->Draw();	
			 				
			PrintHist (histStar, plotName + "60ps", plotTitle, canvasStar);
			
		
		plotName = "HTrackMass";
		plotTitle = "Mass Distribution 20ps";
				
			getCanvas(canvasStar);
			getHist(histFile20ps, plotName, histStar);
			setStyleDistHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetRangeUser(0.01,0.8);
			histStar->GetXaxis()->SetTitle("mass (GeV)");
			histStar->SetStats(0);
			//gStyle->SetOptStat("r");
			histStar->Draw();			
			
			legend = new TLegend(0.1,0.8,0.3,0.9);			
   			legend->AddEntry((TObject*)0, "Time resolution 20ps", "");
   			legend->Draw();	
			 				
	
			 				
			PrintHist (histStar, plotName + "20ps", plotTitle, canvasStar);
			
			
		plotName = "HTrackMass";
		plotTitle = "Mass Distribution 10ps";
				
			getCanvas(canvasStar);
			getHist(histFile10ps, plotName, histStar);
			setStyleDistHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetRangeUser(0.01,0.8);
			histStar->GetXaxis()->SetTitle("mass (GeV)");
			histStar->SetStats(0);
			//gStyle->SetOptStat("r");
			histStar->Draw();
			
			legend = new TLegend(0.1,0.8,0.3,0.9);			
   			legend->AddEntry((TObject*)0, "Time resolution 10ps", "");
   			legend->Draw();	
			 				
	
			 				
			PrintHist (histStar, plotName + "10ps", plotTitle, canvasStar);
			
			
		plotName = "HTrackMass";
		plotTitle = "Mass Distribution 1ps";
				
			getCanvas(canvasStar);
			getHist(histFile1ps, plotName, histStar);
			setStyleDistHist (histStar, plotName, plotTitle, canvasStar);
			histStar->GetXaxis()->SetRangeUser(0.01,0.8);
			histStar->GetXaxis()->SetTitle("mass (GeV)");
			histStar->SetStats(0);
			//gStyle->SetOptStat("r");
			histStar->Draw();
			
			legend = new TLegend(0.1,0.8,0.3,0.9);			
   			legend->AddEntry((TObject*)0, "Time resolution 1ps", "");
   			legend->Draw();	
			 				
	
			 				
			PrintHist (histStar, plotName + "1ps", plotTitle, canvasStar);
			
				
			
//////////////////////////////////////////////////////////////////////////////////////

		
		
			return;
		

//////////////////////////////////////////////////////////////////////////////////////
		}
	
