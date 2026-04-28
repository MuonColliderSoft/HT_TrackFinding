// simulation of track reconstruction. Aug 29 2022
// modified for Muon Collider software integration Apr 14 2025
//
//

#include <iostream> 
#include <fstream> 
#include <iomanip> 
#include <string>
#include <sstream> 
#include <algorithm> 
#include <random>
#include <list>
#include <time.h>

#include "TROOT.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"
//#include "TApplication.h"
#include "TGraph2D.h"
#include <filesystem>

		


using namespace std;

#include "Statistics.cpp"
#include "makeEffGraphFromHists.cpp"
#include "Parameters.h"
#include "CellIDtoLayer.h"
#include "Hit.cpp"
#include "Track.cpp"
#include "HTArray.cpp"
#include "TMath.h"
#include "TTree.h"
#include "TrackReader.h"
#include "BibFileReader.cpp"
#include "Event.cpp"
#include "HTAmapper.h"
#include "CellMap.cpp"


double Pi = 3.14159;


std::string CellIDtoString(int CellID)
{
    const unsigned int system = (unsigned)( CellID & 0x1f );
    const unsigned int side   = (unsigned)( (CellID >> 5)  & 0x3 );
    const unsigned int layer  = (unsigned)( (CellID >> 7)  & 0x3f );
    const unsigned int module = (unsigned)( (CellID >> 13) & 0x7ff );
    const unsigned int sensor = (unsigned)( (CellID >> 24) & 0xff );

    return "Sy" + std::to_string(system)
         + "Si" + std::to_string(side)
         + "La" + std::to_string(layer)
         + "Mo" + std::to_string(module)
         + "Se" + std::to_string(sensor);
}


int TrainingPhase = 0; // 0 = pattern recognition 
                       // 1 = from scratch
                       // 2 = second step after diagonalization
                       // 3 = verification and statistics

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// random generators

//std::default_random_engine generator;
std::mt19937 generator; // Mersenne Twister
std::mt19937 generator_trk; // Mersenne Twister for track generation

std::uniform_real_distribution<double> distribution(0.,1.);
std::normal_distribution<double> gauss(0.,1.);

Parameters par; // class that holds all configurable parameters					   

bool Debug = false;
bool verbose = false;
				   
// The following global parameters are initialized in the main function from Parameters.h

unsigned nEvents; // Number of events to be generated
unsigned nTracks; // Number of tracks per event
long int backGnd; // Number of background events to be simulated


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
 
 
  struct Results {
    	double min;
    	double max;
    	double qmin;
    	double qmax;
    	double cutLow;
    	double cutHigh; 
    	unsigned N = 0;
    	 
    };
    
 
    
    double quantile = 0.90;
    
    void defineCuts(Results &R){
    	double ratioUp = (R.max-R.qmax)/(R.max-R.min)/(1-quantile)*2.;
    		if(ratioUp < 3) R.cutHigh = R.max;
    	
    	double ratioDown = (R.qmin-R.min)/(R.max-R.min)/(1-quantile)*2.;
    		if(ratioDown < 3) R.cutLow = R.min;
    
    };
    




int main(){

	struct X1X2T {
		double x1;
		double x2;
		double t;
    };
    
   
   

	map <int, vector<X1X2T> > SensorToX1X2Tlist;
	map <string, Results> sensorNameToResults;

	// parameter initialization
    nEvents = par.reco_nEvents; // Number of events to be generated
	nTracks = par.reco_nTracks; // Number of tracks per event
	backGnd = par.reco_backGnd; // Number of background events to be simulated
	
	// seeding the random generators 
  	long int randomSeed = par.reco_randomSeed;
  	if(randomSeed)generator.seed(randomSeed);
  	
  	 // Opening the histogram file
  	  	
  	string histFileName = par.reco_histFileName;
	TFile* histFile = new TFile(histFileName.c_str(),"RECREATE");  // histogram file
	
	unsigned nBibHits = 0;
	
	BibFileReader bibRead;
	
	if(backGnd > 0){
	
		// opening BIB background file
	
		cout << "reading BIB... ";
	
		string bibFileName = par.reco_bibFileName;	
		int code = bibRead.readFile(bibFileName);
		cout << "return code = " << code << endl;
		cout << "pool size = " << bibRead.size() << endl;
		if(code == -1 && backGnd > 0) return 0;
	
		unsigned totBIB = bibRead.size();
		cout << totBIB << " BIB events read from file" << endl;	
	
	}
					
	std::poisson_distribution<unsigned> poissDist(nBibHits);
	
	cerr << "------------- init histograms------------------" << endl;
	
		
	histFile->cd();

	
	// MAIN LOOP ON EVENT GENERATION  ///////////////////////////////////////////////
		
	// instantiate track reader from file
	
	TrackReader* newRecoTrack = new TrackReader(par.reco_inputTrackFileName);
	
	long NtotEvents = 0;
	
	cout << "BEGIN EVENT GENERATION " << endl;	

    for(unsigned iEv = 0; iEv != nEvents; ++iEv){
    
		//cout << "********************************************************************" << endl;	
			
	
		// create one event
				
		Event ev(newRecoTrack, nTracks, bibRead, backGnd);// (track reader, nTracks, BIB reader, N bkg hits)
		if((int)ev.trackList.size() < (int)nTracks) break; // end of input file
		
		//if(true) ev.print(cout,1);	
		
		++NtotEvents;
								 
		//if(true) cout << NtotEvents << " events generated" << endl;
	
	
		
		// LOOP ON HITS
	
		unsigned nHits = ev.hitList.size();
		
		//if(verbose) cout << "nHits " << nHits << endl;
		
		for(unsigned iH = 0; iH != nHits; ++iH) {
			
			Hit thisHit = ev.hitList[iH];
			if(!thisHit.isSeed()) continue;// ignore all non-seed hits
			
			//thisHit.print(cout);
			//cout << "***" << endl;
			
			(SensorToX1X2Tlist[thisHit.CellID]).push_back({thisHit.x1,thisHit.x2,thisHit.t});
		
				
		}// end loop on hits	
	
	} // end loop on events ///////////////////////////////////////////////////////////
	  /////////////////////////////////////////////////////////////////////////////////
	
		Statistics statX1, statX2, statT;
									
		// Run ROOT in batch mode (no GUI, safe)
		gROOT->SetBatch(true);		
		// Create reusable histograms
		TH1D h1("h1", "Generic Histogram", 100, 0.0, 1.0);
		h1.SetDirectory(nullptr);  // prevent ROOT ownership issues			
		TH1D h2("h2", "Generic Histogram", 100, 0.0, 1.0);
		h2.SetDirectory(nullptr);  // prevent ROOT ownership issues			
		TH1D ht("ht", "Generic Histogram", 100, 0.0, 1.0);
		ht.SetDirectory(nullptr);  // prevent ROOT ownership issues	
		
			
		// Canvas for drawing
		TCanvas c("c", "canvas", 800, 600);
			
		for (auto& pair : SensorToX1X2Tlist) { // begin loop on sensors
						
				if(pair.second.size() < 1000) continue;					
				string sensorName =  CellIDtoString(pair.first);				
				//std::cout << sensorName << " -> " << pair.second.size();			
					statX1.Clear();
					statX2.Clear();
					statT.Clear();  						
				for (auto i = 0; i != pair.second.size(); ++i){
					statX1.Fill(pair.second[i].x1);
					statX2.Fill(pair.second[i].x2);
					statT.Fill(pair.second[i].t);  		
				}
				//cout << " mx1 = " << statX1.GetMean() << " sx1 = " << statX1.GetSigma() 
				//		<< " min = " <<  statX1.GetMin() << " max = " <<  statX1.GetMax();
				//cout << " mx2 = " << statX2.GetMean() << " sx2 = " << statX2.GetSigma() 
				//		<< " min = " <<  statX2.GetMin() << " max = " <<  statX2.GetMax(); 
				//cout << " mT = " << statT.GetMean() << " sT = " << statT.GetSigma() 
				//		<< " min = " <<  statT.GetMin() << " max = " <<  statT.GetMax() << endl; 												
		///////////////////////////////////////////////////////////////////////////////////		    			
		///////////////////////////////////////////////////////////////////////////////////		
		// GENERATING HISTOGRAMS 
		///////////////////////////////////////////////////////////////////////////////////		    			
		///////////////////////////////////////////////////////////////////////////////////	
		
				// Make sure output folder exists (creates it if missing)
				std::filesystem::create_directories("hists");	
		
				// --- Define parameters for this iteration ---
				
				double xmin = statX1.GetMin();
				double xmax = statX1.GetMax();
		
				std::string title    = sensorName + "X1";
				std::string xlabel   = "X1";
				std::string filename = "hists/" + sensorName + "_X1.pdf";
		
				// --- Reset histogram contents ---
				h1.Reset();
		
				// --- Redefine binning (this sets X range properly) ---
				h1.SetBins(100, xmin, xmax);
		
				// --- Set titles ---
				h1.SetTitle(title.c_str());
				h1.GetXaxis()->SetTitle(xlabel.c_str());
				h1.GetYaxis()->SetTitle("Entries");
		
				// --- Fill histogram ---
				for (auto i = 0; i != pair.second.size(); ++i){
					h1.Fill(pair.second[i].x1);		
				}
				
				histFile->cd();

				std::string rootName = sensorName + "_X1";
				TH1D* hsave = (TH1D*)h1.Clone(rootName.c_str());
				hsave->Write();
				delete hsave;
		
				// --- Draw and save ---
				c.cd();
				//c.SetLogy();
				h1.Draw();
				c.Print(filename.c_str());
		
				//std::cout << "Saved " << filename << std::endl;
				
				double qmin, qmax;
				statX1.GetQuantile(quantile, qmin,qmax);
				//cout << "x1max: " << xmax << " x1min: " << xmin << " qmin: " << qmin << " qmax: " << qmax << endl;		
				double cutLow = qmin;
				double cutHigh = qmax;
				
				//store results in map
				sensorNameToResults[sensorName + "_X1"].min = xmin;
				sensorNameToResults[sensorName + "_X1"].max = xmax;
				sensorNameToResults[sensorName + "_X1"].qmin = qmin;
				sensorNameToResults[sensorName + "_X1"].qmax = qmax;		
				sensorNameToResults[sensorName + "_X1"].cutLow = cutLow;						
				sensorNameToResults[sensorName + "_X1"].cutHigh = cutHigh;
					
		///////////////////////////////////////////////////////////////////////////////////					
			
				// --- Define parameters for this iteration ---
			{
				 xmin = statX2.GetMin();
				 xmax = statX2.GetMax();
		
				 title    = sensorName + "X2";
				std::string xlabel   = "X2";
				std::string filename = "hists/" + sensorName + "_X2.pdf";
		
				// --- Reset histogram contents ---
				h2.Reset();
		
				// --- Redefine binning (this sets X range properly) ---
				h2.SetBins(100, xmin, xmax);
		
				// --- Set titles ---
				h2.SetTitle(title.c_str());
				h2.GetXaxis()->SetTitle(xlabel.c_str());
				h2.GetYaxis()->SetTitle("Entries");
		
				// --- Fill histogram ---
				for (auto i = 0; i != pair.second.size(); ++i){
					h2.Fill(pair.second[i].x2);  		
				}
				
				histFile->cd();

				std::string rootName = sensorName + "_X2";
				TH1D* hsave = (TH1D*)h2.Clone(rootName.c_str());
				hsave->Write();
				delete hsave;
		
				// --- Draw and save ---
				c.cd();
				//c.SetLogy();
				h2.Draw();
				c.Print(filename.c_str());
		
				//std::cout << "Saved " << filename << std::endl;
				
				statX2.GetQuantile(quantile, qmin,qmax);
				//cout << "x2max: " << xmax << " x2min: " << xmin << " qmin: " << qmin << " qmax: " << qmax << endl;
				 cutLow = qmin;
				 cutHigh = qmax;
				
				//store results in map
				sensorNameToResults[sensorName + "_X2"].min = xmin;
				sensorNameToResults[sensorName + "_X2"].max = xmax;
				sensorNameToResults[sensorName + "_X2"].qmin = qmin;
				sensorNameToResults[sensorName + "_X2"].qmax = qmax;		
				sensorNameToResults[sensorName + "_X2"].cutLow = cutLow;						
				sensorNameToResults[sensorName + "_X2"].cutHigh = cutHigh;
					
			}	    			
		///////////////////////////////////////////////////////////////////////////////////					
			{			
				// --- Define parameters for this iteration ---
				double xmin = statT.GetMin();
				double xmax = statT.GetMax();
		
				std::string title    = sensorName + "T";
				std::string xlabel   = "T";
				std::string filename = "hists/" + sensorName + "_T.pdf";
		
				// --- Reset histogram contents ---
				ht.Reset();
		
				// --- Redefine binning (this sets X range properly) ---
				ht.SetBins(100, xmin, xmax);
		
				// --- Set titles ---
				ht.SetTitle(title.c_str());
				ht.GetXaxis()->SetTitle(xlabel.c_str());
				ht.GetYaxis()->SetTitle("Entries");
		
				// --- Fill histogram ---
				for (auto i = 0; i != pair.second.size(); ++i){
					ht.Fill(pair.second[i].t);		
				}
				
				histFile->cd();

				std::string rootName = sensorName + "_T";
				TH1D* hsave = (TH1D*)ht.Clone(rootName.c_str());
				hsave->Write();
				delete hsave;
		
				// --- Draw and save ---
				c.cd();
				//c.SetLogy();
				ht.Draw();
				c.Print(filename.c_str());
		
				//std::cout << "Saved " << filename << std::endl;
				
				statT.GetQuantile(quantile, qmin,qmax);
				//cout << "tmax: " << xmax << " tmin: " << xmin << " qmin: " << qmin << " qmax: " << qmax << endl;
				 cutLow = qmin;
				 cutHigh = qmax;
				
				//store results in map
				sensorNameToResults[sensorName + "_T"].min = xmin;
				sensorNameToResults[sensorName + "_T"].max = xmax;
				sensorNameToResults[sensorName + "_T"].qmin = qmin;
				sensorNameToResults[sensorName + "_T"].qmax = qmax;		
				sensorNameToResults[sensorName + "_T"].cutLow = cutLow;						
				sensorNameToResults[sensorName + "_T"].cutHigh = cutHigh;
		
			}	    			
		///////////////////////////////////////////////////////////////////////////////////					
							
									
						 
	} // end loop on sensors
	
	
	
 	cout << NtotEvents << " events generated" << endl;
  
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

	histFile->Write();
	histFile->Close();
	delete histFile;
	

// ------------------------------------------------------------
// Print formatted results table with row/column separators
// ------------------------------------------------------------

cout << std::fixed << std::setprecision(1);

const int Wname = 28;
const int Wnum  = 10;

auto line = [&]() {
    cout << "+"
         << std::string(Wname, '-')
         << "+"
         << std::string(Wnum, '-')
         << "+"
         << std::string(Wnum, '-')
         << "+"
         << std::string(Wnum, '-')
         << "+"
         << std::string(Wnum, '-')
         << "+"
         << std::string(Wnum, '-')
         << "+"
         << std::string(Wnum, '-')
         << "+"
         << "\n";
};

cout << "\n";
line();

cout << "|"
     << std::setw(Wname) << "Sensor"
     << "|"
     << std::setw(Wnum) << "Min"
     << "|"
     << std::setw(Wnum) << "Max"
     << "|"
     << std::setw(Wnum) << "QMin"
     << "|"
     << std::setw(Wnum) << "QMax"
     << "|"
     << std::setw(Wnum) << "CutLow"
     << "|"
     << std::setw(Wnum) << "CutHigh"
     << "|"
     << "\n";

line();

// Rows
for (auto& pair : sensorNameToResults) {

	defineCuts(pair.second); // define cuts according to quantiles and min,max

    cout << "|"
         << std::setw(Wname) << pair.first
         << "|"
         << std::setw(Wnum) << pair.second.min
         << "|"
         << std::setw(Wnum) << pair.second.max
         << "|"
         << std::setw(Wnum) << pair.second.qmin
         << "|"
         << std::setw(Wnum) << pair.second.qmax
         << "|"
         << std::setw(Wnum) << pair.second.cutLow
         << "|"
         << std::setw(Wnum) << pair.second.cutHigh
         << "|"
         << "\n";

    line();
}

// READ BACKGROUND HITS

unsigned NTotBkg;

struct counts {
	unsigned NX1 = 0;
	unsigned NX2 = 0;
	unsigned NT  = 0;
	unsigned N   = 0;
	unsigned Nin = 0;
};

map<string, counts>  sensorNameToAcceptedCount;

{	
	cout << "Reading background..." << endl;
	
	string BibFileName = "ntu_bkg_hits_phim22-52-1evt.root";
	TrackReader* bibReader = new TrackReader(BibFileName);
	Track thisTrack;
	
	bool eof = bibReader->read(thisTrack);
		if(eof) {
			cerr << "BibFileReader failed" << endl;
			
		}		
		// loop on hits of this track
		NTotBkg = (int)thisTrack.hitList.size();
		cout << NTotBkg << " reading background hits" << endl;
		cout << "Please wait...." << endl << std::flush;
							
		for(int iH = 0; iH != (int)thisTrack.hitList.size(); ++iH){
			Hit thisHit = thisTrack.hitList[iH];
			
			// get the sensor name for this hit
			string sensorName =  CellIDtoString(thisHit.CellID);
			
			bool X1cut = false;
			bool X2cut = false;
			bool Tcut  = false;
						
											
			auto it1 = sensorNameToResults.find(sensorName + "_X1");		
			if (it1 != sensorNameToResults.end()){
				++sensorNameToAcceptedCount[sensorName].Nin; // counts per sensor		
				//cout << "X1 FOUND ";
				Results& R = it1->second;		
				X1cut = (thisHit.x1 >= R.cutLow && thisHit.x1 <= R.cutHigh);
				//cout << X1cut << endl;
				if(X1cut) ++sensorNameToAcceptedCount[sensorName].NX1;
			}
			auto it2 = sensorNameToResults.find(sensorName + "_X2");													
			if (it2 != sensorNameToResults.end()){
				//cout << "X2 FOUND ";
				Results& R = it2->second;		
				X2cut = (thisHit.x2 >= R.cutLow && thisHit.x2  <= R.cutHigh);
				//cout << X2cut << endl;
				if(X2cut) ++sensorNameToAcceptedCount[sensorName].NX2;
			}		
			auto it3 = sensorNameToResults.find(sensorName + "_T");													
			if (it3 != sensorNameToResults.end()){
				//cout << "T FOUND ";
				Results& R = it3->second;		
				 Tcut = (thisHit.t >= R.cutLow && thisHit.t <= R.cutHigh);
				 //cout << Tcut << endl;
				 if(Tcut) ++sensorNameToAcceptedCount[sensorName].NT;		
				 if (X1cut && X2cut && Tcut) {
				 ++sensorNameToAcceptedCount[sensorName].N;
				}
			}	
			
		} // end loop on hits of this track
						
}		cout << "Finished reading background hits" << endl;
		

// ------------------------------------------------------------
// FINAL TABLE: accepted background hits per sensor
// ------------------------------------------------------------

cout << std::fixed << std::setprecision(0);

const int Wname2 = 28;
const int Wcnt   = 10;

auto line2 = [&]() {
    cout << "+"
         << std::string(Wname2, '-')
         << "+"
         << std::string(Wcnt, '-')
         << "+"
         << std::string(Wcnt, '-')
         << "+"
         << std::string(Wcnt, '-')
         << "+"
         << std::string(Wcnt, '-')
         << "+"
         << std::string(Wcnt, '-')
         << "+\n";
};

cout << "\n";
cout << "BACKGROUND HITS PASSING CUTS\n";
cout << "(individual coordinates and all three together)\n";

line2();

cout << "|"
     << std::setw(Wname2) << "Sensor"
     << "|"
     << std::setw(Wcnt)   << "Nin  "
     << "|"
     << std::setw(Wcnt)   << "NX1"
     << "|"
     << std::setw(Wcnt)   << "NX2"
     << "|"
     << std::setw(Wcnt)   << "NT"
     << "|"
     << std::setw(Wcnt)   << "NAll"
     << "|\n";

line2();

unsigned totalNin = 0;
unsigned totalNX1 = 0;
unsigned totalNX2 = 0;
unsigned totalNT  = 0;
unsigned totalN   = 0;
unsigned nSensors = 0;

for (const auto& pair : sensorNameToAcceptedCount) {

    const std::string& sensorName = pair.first;
    const counts& C = pair.second;

    cout << "|"
         << std::setw(Wname2) << sensorName
         << "|"
         << std::setw(Wcnt)   << C.Nin
         << "|"
         << std::setw(Wcnt)   << C.NX1
         << "|"
         << std::setw(Wcnt)   << C.NX2
         << "|"
         << std::setw(Wcnt)   << C.NT
         << "|"
         << std::setw(Wcnt)   << C.N
         << "|\n";

    line2();
    
	totalNin += C.Nin;
    totalNX1 += C.NX1;
    totalNX2 += C.NX2;
    totalNT  += C.NT;
    totalN   += C.N;
    ++nSensors;
}

// Totals row
cout << "|"
     << std::setw(Wname2) << "TOTAL"
     << "|"
     << std::setw(Wcnt)   << totalNin
     << "|"
     << std::setw(Wcnt)   << totalNX1
     << "|"
     << std::setw(Wcnt)   << totalNX2
     << "|"
     << std::setw(Wcnt)   << totalNT
     << "|"
     << std::setw(Wcnt)   << totalN
     << "|\n";

line2();

cout << nSensors << " sensors listed.\n";


cout << "PREDICTIONS" << endl;
cout << endl;


cout << std::scientific << std::setprecision(6);
for (const auto& pair : sensorNameToAcceptedCount) {
	const std::string& sensorName = pair.first;
    const counts& C = pair.second;
    
    unsigned NTot = C.Nin;
    
    double pX1 = (double)C.NX1/(double)NTot;
    double pX2 = (double)C.NX2/(double)NTot;
    double pT = (double)C.NT/(double)NTot;
    double PtotExp = pX1*pX2*pT;
    unsigned Nexp = (double)NTot*PtotExp;
    cout << sensorName << "  pX1 " << pX1 << " pX2 " << pX2 << " pT " << pT << " Ptot " << PtotExp << " Nexp: " << Nexp << " actual: " << C.N << endl;

}


	return 0;	
	
		
	
} // end main



