/*

#include <iostream> 
#include <fstream> 
#include <iomanip> 
#include <string>
#include <sstream> 
#include <algorithm> 
#include <random>
#include <list>
#include <time.h>

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph2D.h"

using namespace std;

#include "Statistics.cpp"
#include "Geometry.cpp".q
.
Geometry g;


*/


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////





int PlotTracks(){

	TGraph2D* g2;
	TCanvas* c;


	const int n_max_points = 100;
  	double xp[n_max_points];
  	double yp[n_max_points];
  	double zp[n_max_points];

  
  

	// Open data file
	
	string dataFileName = "PlotData.txt";
	cout << "Opening Data File: " << dataFileName << " ..." << endl;
	ifstream infile;
	infile.open(dataFileName);
	if(infile) cout << " OK " << endl; 
	else {
		cout << "Data File not found, abort." << endl;
		return 0;
	}
	
	// begin loop to read file to the end
	
	string tag;
	unsigned evNum;
	unsigned nCand;
	
	for(;infile;){
	
	
		infile >> tag >> evNum >> nCand;
		if(!infile) return 0;

		//cout << "read " << tag  << " " << evNum << " " << nCand << endl;

		for (int iCand = 0; iCand != nCand; ++ iCand){
			int iPhi, iEta, iInvPt, nLay;
			infile >> iPhi >> iEta >> iInvPt >> nLay;
			//cout << "read " << iPhi  << " " << iEta << " " << iInvPt << " " << nLay << endl;
			cout << endl;
			for(int iLay = 0; iLay != nLay; ++iLay){// begin loop on layers
				int nHits;
				infile >> nHits; //cout << "read " << nHits << endl;
				for(int ih = 0; ih != nHits; ++ih){// begin loop on hits
					char hitType; 
					int layerInd;
					double X,Y,Z,T;
					int trackInd;
					infile >> hitType >> layerInd >> X >> Y >> Z >> T >> trackInd; 
					cout << hitType << " " << layerInd << " "  << X << " "  << Y << " "  << Z << " "  << T << " "  << trackInd << endl; 
					xp[iLay] = X;
					yp[iLay] = Y;
					zp[iLay] = Z;
				}// end loop on hits
				
			}// end loop on layers
			
			
			cout << endl << tag  << " " << evNum << " Candidate: " << iCand << " cell: [" << iPhi << "," << iEta << "," << iInvPt << "]" << endl;
			
			c = new TCanvas("c","Track plot",0,0,800,1250);
			
			g2 = new TGraph2D(nLay, xp, yp, zp);
			
			g2->SetMarkerStyle(20);
  			g2->SetMarkerSize(1);
  
  			TH3I *Track_Plot = new TH3I("Track_Plot","Track_Plot",1,0.,+1600.,1,0.,+1600.,1,0.,+2500.);
  			Track_Plot->SetStats(false);
  
 			Track_Plot->SetTitle("; X [mm];Y [mm];Z [mm]");
  
  			Track_Plot->Draw();
  
  			g2->Draw("PSAME");
			
			
			
			//return 0;
			
			
			//c->ForceUpdate();
			//c->Flush();
			//c->Draw();
			gPad->Update();
			gPad->WaitPrimitive();
			
		/*	cout << "next [y/n]?";
			string cc;
			cin >> cc;
			if(cc == "n") return 0;	
		*/	
			delete c;
			delete g2;
			delete Track_Plot;
		
			
		}// end loop on candidates
	}// end loop on events
	
	return 0;
	

	
} // end main



