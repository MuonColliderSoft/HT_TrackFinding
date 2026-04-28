// simulation of track reconstruction. Aug 29 2022
// modified for Muon Collider software integration Apr 14 2025
//
//==============================================================================
// Reconstruction.cpp
//
// Purpose:
//   Muon Collider track-reconstruction study using a Hough Transform Array
//   (HTA). Reads simulated signal tracks, optionally overlays background hits,
//   performs pattern recognition, fits track candidates, and stores
//   performance histograms.
//
// Main Steps:
//   1. Initialize HTA and load training data.
//   2. Read events with signal + optional BIB background.
//   3. Fill HTA with seed hits.
//   4. Find candidate cells.
//   5. Fit candidate tracks.
//   6. Apply chi2 quality cuts.
//   7. Remove duplicate candidates.
//   8. Write ROOT histograms.
//
// Outputs:
//   - ROOT histogram file
//   - Candidate / fit statistics
//   - Signal and background hit distributions
//   - Reconstruction performance metrics
//
// Example Build:
//   clang++ -std=c++17 -O2 -Wall -Wextra Reconstruction.cpp -o Reconstruction \
//       `root-config --cflags --libs`
//
// Example Run:
//   ./Reconstruction
//
// Requirements:
//   ROOT, C++17, Muon Collider support classes.
//
//==============================================================================

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
#include "TGraph2D.h"


using namespace std;

//#include "Statistics.cpp"
//#include "makeEffGraphFromHists.cpp"
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

bool Debug = false;
bool Special = false;
bool verbose = false;
				   
// The following global parameters are initialized in the main function from Parameters.h

unsigned nEvents; // Number of events to be generated
unsigned nTracks; // Number of tracks per event
long int backGnd; // Number of background events to be simulated
int fillMode; // optimization mode for HTA fill
bool PlotTracks; // create a file with data for 3-D plots of candidates	


int NFITS = 0;
int TOTNFITS = 0;

	
Parameters par; // global class that holds all configurable parameters					   


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

int main(){

	
	
	unsigned HTA_NphiBins = par.HTA_NphiBins;
	unsigned HTA_NetaBins = par.HTA_NetaBins;
	unsigned HTA_NinvptBins = par.HTA_NinvptBins;

	// parameter initialization
    nEvents = par.reco_nEvents; // Number of events to be generated
	nTracks = par.reco_nTracks; // Number of tracks per event
	backGnd = par.reco_backGnd; // Number of background events to be simulated
	fillMode = par.reco_fillMode; // optimization mode for HTA fill
	PlotTracks = par.reco_PlotTracks; // create a file with data for 3-D plots of candidates
	
	// seeding the random generators 
  	long int randomSeed = par.reco_randomSeed;
  	if(randomSeed)generator.seed(randomSeed);
  	
  	
  	 // Opening the histogram file
  	  	
  	string histFileName = par.reco_histFileName;
	TFile* histFile = new TFile(histFileName.c_str(),"RECREATE");  // histogram file
	

	///////////////////////////////////////////////////////////////////////////	
	// Begin Instantiate Hough Transform Array ////////////////////////////////
	
	cout << "Instantiating Hough Transform Array..." << endl;
	if(verbose)cout << "verbose..." << endl;
	cout << "***"<< endl;
		
	HTArray *HTA = new HTArray();
		
	HTA->initHists();
	HTA->print(cout);
	
	// Check sanity of array dimensions
	// see CellMap.cpp
	
	IndexCodec::checkDimensions();
			
	// Open training data file for HTA
	
	string dataFileName = par.reco_dataFileName;
	cout << "Opening Training Data File: " << dataFileName << " ..." << endl;
	ifstream infile;
	infile.open(dataFileName);
	
	if(infile){ // data file found
		int retcode = HTA->read(infile);
		cout << "HTA.read retcode: " << retcode << endl;
					
		if(retcode < 0) {
			cout << "datafile does not match current configuration" << endl;
			return retcode; 
		}
		else cout << "HTA initialized correctly" << endl;		
	} 	
	else {
		cout << "Data File not found, abort." << endl;
		return 0;
	}
	// End of Hough Transform Array instantiation ////////////////////////////////	
	//////////////////////////////////////////////////////////////////////////////
	
	ofstream outPlotFile; // file for plotting candidates
		
	if(PlotTracks) {		
		// Open plot tracks data file	
		string plotDataFileName = par.reco_plotDataFileName;
		cout << "Opening Plot Data File: " << plotDataFileName << " ..." << endl;	
		outPlotFile.open(plotDataFileName);	
		if(outPlotFile){ // plot data file opened successfully
			cout << plotDataFileName << " opened successfully "<< endl;					
		} 	
		else {
			cout << "Error opening " << plotDataFileName << endl;
			return 0;
		}	
	}	
	/////////////////////////////////////////////////////////////////////////////////////////
	
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
	
	
	
/*	TH1D HnHitsInThisCell("HnHitsInThisCell","HnHitsInThisCell", 13, -0.5,+12.5);
	TH1D HnBestCellHits("HnBestCellHits","HnBestCellHits", 13, -0.5,+12.5);
	TH1D HCellStat("HCellStat","HCellStat",21, -0.5, 20.5);
*/	

	cerr << "------------- init histograms------------------" << endl;
	
		
	histFile->cd();
	
	
				
	TH1D HnCandidates("HnCandidates","HnCandidates", 501, -0.5,+500.5);
	TH1D HnFits("HnFits","HnFits", 1001,-0.5, 1000.5);
	
	TH1D HSighitx("HSighitx","HSighitx",1000.,-1500.,+1500);
	TH1D HSighity("HSighity","HSighity",1000.,-1500.,+1500);
	TH1D HSighitz("HSighitz","HSighitz",1000.,-1500.,+1500);
	TH1D HSighitt("HSighitt","HSighitt",1000.,-1500.,+1500);
	TH2I H2Sighitxy("H2Sighitxy","H2Sighitxy",1000,-1500., +1500.,1000,-1500,+1500);		
	TH2I H2SighitzR("H2SighitzR","H2SighitzR",1000,-1500., +1500.,1000,-0.,+1500);
	
	TH1D HBkghitx("HBkghitx","HBkghitx",1000.,-1500.,+1500);
	TH1D HBkghity("HBkghity","HBkghity",1000.,-1500.,+1500);
	TH1D HBkghitz("HBkghitz","HBkghitz",1000.,-1500.,+1500);
	TH1D HBkghitt("HBkghitt","HBkghitt",1000.,-1500.,+1500);
	TH2I H2Bkghitxy("H2Bkghitxy","H2Bkghitxy",1000,-1500., +1500.,1000,-1500,+1500);		
	TH2I H2BkghitzR("H2BkghitzR","H2BkghitzR",1000,-1500., +1500.,1000,-0.,+1500);
	
	// CHI SQUARES PER N HITS ////////////////////////////////////
	
	TH1D HFitLayers("HFitLayers","HFitLayers",21, -0.5, 20.5);
	TH1D HFitAllChi2("HFitAllChi2","HFitAllChi2",5000, 0., 5000.);
	TH1D HFitGoodChi2("HFitGoodChi2","HFitGoodChi2",5000, 0., 5000.);
	TH1D HFitRedChi2("HFitRedChi2","HFitRedChi2",5000, 0., 5000.);
/*	
	TH1D HFit5Chi2("HFit5Chi2","HFit5Chi2",5000, 0., 5000.);
	TH1D HFit6Chi2("HFit6Chi2","HFit6Chi2",5000, 0., 5000.);
	TH1D HFit7Chi2("HFit7Chi2","HFit7Chi2",5000, 0., 5000.);
	TH1D HFit8Chi2("HFit8Chi2","HFit8Chi2",5000, 0., 5000.);	
	TH1D HFit9Chi2("HFit9Chi2","HFit9Chi2",5000, 0., 5000.);
	TH1D HFit10Chi2("HFit10Chi2","HFit10Chi2",5000, 0., 5000.);
	TH1D HFit11Chi2("HFit11Chi2","HFit11Chi2",5000, 0., 5000.);
	TH1D HFit12Chi2("HFit12Chi2","HFit12Chi2",5000, 0., 5000.);
	TH1D HFit13Chi2("HFit13Chi2","HFit13Chi2",5000, 0., 5000.);
	TH1D HFit14Chi2("HFit14Chi2","HFit14Chi2",5000, 0., 5000.);
	TH1D HFit15Chi2("HFit15Chi2","HFit15Chi2",5000, 0., 5000.);
	TH1D HFit16Chi2("HFit16Chi2","HFit15Chi2",5000, 0., 5000.);
*/	
	

	
	TH1D HnFoundTracks("HnFoundTracks","HnFoundTracks",3, -0.5, 2.5);
	
	

	
	// MAIN LOOP ON EVENT GENERATION  ///////////////////////////////////////////////
	
	
	struct FitTrack {
	
        		double chi2;
        		double phi;
        		double eta;
        		double invPt; 
        		double z0;
        		double t0;
        		double beta; 
        		unsigned nLayers;
    };
  	
	
	
	double candidateRate = 0.;
	unsigned nEventsWithCandidates = 0;
	unsigned nEventsWithTracks = 0;
	
	// instantiate track reader from file
	
	TrackReader* newRecoTrack = new TrackReader(par.reco_inputTrackFileName);
	
	
	//int argc; 
	//char **argv;

//	TApplication theApp("app",&argc, argv);
	
	//TGraph2D* g2;
  	//TCanvas* c;
 
	
	cout << "BEGIN EVENT GENERATION " << endl;
	
	time_t t0     = time(NULL);  // start of loop (total time)
	time_t t_prev = t0;          // previous iteration time
	
	
	// begin loop on events
	
	for(unsigned iEv = 0; iEv != nEvents; ++iEv){
	
		cout << "**********************************************************************************************" << endl;
		cout << "**********************************************************************************************" << endl;
		cout << "**********************************************************************************************" << endl;
	
		// current time
		time_t my_time = time(NULL);
	
		// print human-readable time
		printf("%s", ctime(&my_time));
	
		// elapsed since previous iteration
		double delta = difftime(my_time, t_prev);
	
		// elapsed since beginning
		double total = difftime(my_time, t0);
	
		cout << iEv << "/" << nEvents << " processed events"
			 << "   (elapsed time: " << delta << " s"
			 << ", total: " << total << " s)" << endl;
	
		// update for next iteration
		t_prev = my_time;

		NFITS = 0;
	
		// create one event
		
		Event ev(newRecoTrack, nTracks, bibRead, backGnd);// (track reader, nTracks, BIB reader, N bkg hits)
		if((int)ev.trackList.size() < (int)nTracks) break; // end of input file
		
		if(true) ev.print(cout,1);	
								 
		if(verbose) cout << "event generated" << endl;
		vector <FitTrack> foundTracks; // tracks found in this event (now empty)
		

		HTA->reset(); // clear all hits from HT Array
		
		if(verbose) cout << "HTA reset" << endl;
		

		
		// LOOP ON HITS
	
		unsigned nHits = ev.hitList.size();
		
		if(verbose) cout << "nHits " << nHits << endl;
		
		for(unsigned iH = 0; iH != nHits; ++iH) {
			
			Hit thisHit = ev.hitList[iH];
			if(!thisHit.isSeed()) continue;// ignore all non-seed hits
			
				double x = thisHit.x;
				double y = thisHit.y;
				double z = thisHit.z;
				double t = thisHit.t;
				
				double R = sqrt(x*x + y*y);
				
			if(thisHit.trackInd == 0) { //this id background			
				HBkghitx.Fill(x);
				HBkghity.Fill(y);
				HBkghitz.Fill(z);
				HBkghitt.Fill(t);				
				H2Bkghitxy.Fill(x,y);		
				H2BkghitzR.Fill(z,R);			
			}
			else{					// this is signal
				HSighitx.Fill(x);
				HSighity.Fill(y);
				HSighitz.Fill(z);
				HSighitt.Fill(t);				
				H2Sighitxy.Fill(x,y);							
				H2SighitzR.Fill(z,R);				
			}
			
			
		
			//////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////			
			//////////////////////////////////////////////////////////////////////////////
			/// FILL HTA ARRAY FOR PATTERN REGOGNITION ///////////////////////////////////
		
			
			if(thisHit.isSeed())HTA->fill(thisHit, par.reco_fillMode);// use only seed hits 
						
		}// end loop on hits
		
		unsigned phi_b, eta_b, pt_b;// best cell
		
		unsigned nLayersHit = HTA->getBestCell(phi_b, eta_b, pt_b);// num of layers hit in best cell
		
		cout << endl << "bestCell: [" << phi_b << "," << eta_b << "," << pt_b << "] nLayersHit: " << nLayersHit << endl;	
		HTA->ArrElem[phi_b][eta_b][pt_b].printHits(cout);
			
		
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
				
		// Deal with candidates
		
		
				
			// Deal with candidates
							
			unsigned nCan = HTA->getCellCandidates();
			
			cout << "N Candidates: " << nCan << endl;
			
			HnCandidates.Fill(nCan);			
			if(nCan) ++nEventsWithCandidates;
		
			bool foundTrack = false;
			
			// Loop on candidates 
			
			for(unsigned iCan = 0; iCan != nCan; ++iCan){
			
				// TrackFit	
			
				HTArray::Pars p = HTA->cellCandidateList[iCan];	
					
				if(verbose) cout << "fit candidate " << iCan << " [" 
					<< p.iPhi << "][" << p.iEta << "][" << p.iInvpt << "]" << endl;					
				double chi2, phi, eta, invPt, z0, t0, beta;
				vector <long int> goodFitHitList;
				unsigned nHitsFit;
						
				int retCodeFit = HTA->ArrElem[p.iPhi][p.iEta][p.iInvpt].
						fitCandidate(&(HTA->allHits),chi2, phi, eta, invPt, z0, t0, beta, nHitsFit, goodFitHitList);
				
				//HnFits.Fill(NFITS);
					
				if(retCodeFit != 0){
					cout << "FIT FAILED with retCodeFit = " << retCodeFit << endl;			
					continue;		
				}	
								
				unsigned dof = 3*nHitsFit - 5; // Degrees of freedom of the fit			
				cout << endl << "retCodeFit: " << retCodeFit << " chi2 = "<< chi2 << " d.o.f. = " << dof << endl;								
				double redChi2 = chi2/double(dof);// reduced chi2						
				HFitLayers.Fill(nHitsFit);
				HFitAllChi2.Fill(chi2);		
				HFitRedChi2.Fill(redChi2);		
				
				Track thisTrack = ev.trackList[0];
				cout << endl << "Generated track: phi = " << thisTrack.phi << " eta = " << thisTrack.eta << " invPt = " << thisTrack.invPt 
					<< " z0 = " << thisTrack.z0 << " t0 = " << thisTrack.t0 << " beta = " << thisTrack.beta << endl;
										
				cout << "Fit result:      phi = " << phi << " eta = " << eta << " invPt = " << invPt 
					<< " z0 = " << z0 << " t0 = " << t0 << " beta = " << beta << endl;
		
				cout << "Number of fits: " << NFITS << endl;
				cout <<"Total Number of Fits = " << TOTNFITS << endl;								
				
								
				if( (retCodeFit == 0) && (redChi2 <= par.reco_chi2Cut) ) {//this is a good fit
				
					cout << "This is a good fit. Found a good track." << endl;
																	
					foundTrack = true;
					HFitGoodChi2.Fill(redChi2);	
					
					// trick to eliminate duplicates							
					
					int i = phi_to_xi(phi);
					int j = eta_to_xj(eta);
					int k = invPt_to_xk(invPt);
						
						if(i >= HTA_NphiBins) i = HTA_NphiBins - 1;
						if(j >= HTA_NetaBins) j = HTA_NetaBins - 1;
						if(k >= HTA_NinvptBins) k = HTA_NinvptBins - 1;
						
						if(i < 0) i = 0;
						if(j < 0) j = 0;
						if(k < 0) k = 0;
		
					bool done =  HTA->ArrElem[i][j][k].thisCellDone;
					
					if(!done){ // this is executed only for one copy of duplicates
					
						HTA->ArrElem[i][j][k].thisCellDone = true;
												
						HFitLayers.Fill(nHitsFit);
											
						FitTrack ft;
							ft.chi2 = chi2;
							ft.phi = phi;
							ft.eta = eta;
							ft.invPt = invPt;
							ft.z0 = z0;
							ft.t0 = t0;
							ft.beta = beta;
							ft.nLayers = nHitsFit;
																				
						foundTracks.push_back(ft);

			   
							   
			
						if(verbose) cout << " good fit"	 << endl;					
										
						int DoF = 3*nHitsFit - 5;			
						if(verbose) cout << " DoF " << DoF << " chi2 " << chi2 << endl;
						if(verbose) cout << " result  pt:" << 1./invPt << " eta: " << eta << " phi: " << phi << " z0: " << z0 << " t0:" << t0 << endl;					
						if(ev.trackList.size()) 
							if(verbose) cout << "  track  pt:" << 1./ev.trackList[0].invPt << " eta: " << ev.trackList[0].eta 
							<< " phi: " << ev.trackList[0].phi << " z0: " << ev.trackList[0].z0 << " t0:" << ev.trackList[0].t0 << endl;
					}// end trick to eliminate duplicates
											
				}// end this is a good fit (retcode and chi2)
												
							
			} // end loop on candidates
		
		cout << endl << foundTracks.size() << " Tracks found for this event" << endl;
		HnFoundTracks.Fill(foundTracks.size());
		
	} // end loop on events ///////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////
	
		
	cout << "Writing histogram file " << histFileName << "..." << endl;
	histFile->Write(); // write histogram file
	
	return 0;	
	
	
		
	
} // end main



