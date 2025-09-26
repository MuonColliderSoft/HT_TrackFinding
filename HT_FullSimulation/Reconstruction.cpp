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
#include "TApplication.h"
#include "TGraph2D.h"


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


//#include "Gauss.h"


//#include "GlobalConstants.cpp"

//#include "Parameters.h"
Parameters par; // class that holds all configurable parameters					   

//#include "DetectorGeometry.cpp"
//DetectorGeometry dg;

//#include "TrackGeometry.cpp"
//TrackGeometry tg("gen"); // track parameters for event generation

bool Debug = false;
bool verbose = par.reco_verbose;
bool Special = false;


				   
// The following global parameters are initialized in the main function

unsigned nEvents; // Number of events to be generated
unsigned nTracks; // Number of tracks per event
long int backGnd; // Number of background events to be simulated
//int fillMode; // optimization mode for HTA fill
bool PlotTracks; // create a file with data for 3-D plots of candidates		


				   

int main(){

	// parameter initialization
    nEvents = par.reco_nEvents; // Number of events to be generated
	nTracks = par.reco_nTracks; // Number of tracks per event
	backGnd = par.reco_backGnd; // Number of background events to be simulated
	int fillMode = par.reco_fillMode; // optimization mode for HTA fill
	PlotTracks = par.reco_PlotTracks; // create a file with data for 3-D plots of candidates
	
	// seeding the random generators 
  	long int randomSeed = par.reco_randomSeed;
  	if(randomSeed)generator.seed(randomSeed);
  	
  	string histFileName = par.reco_histFileName;
	TFile* histFile = new TFile(histFileName.c_str(),"RECREATE");  // histogram file
  
	// Instantiate Hough Transform Array
	
	cout << "Instantiating Hough Transform Array..." << endl;
	
	
	HTArray *HTA = new HTArray();
	
	HTA->initHists();
	HTA->print(cout);
	
	// Check sanity of array dimensions
	// see CellMap.cpp
	
	IndexCodec::checkDimensions();
			
	// Open data file
	
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
	} 
	
	else {
		cout << "Data File not found, abort." << endl;
		return 0;
	}
		
	
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
	
	
/*			
	TH1D HnCandidates("HnCandidates","HnCandidates", 6, -0.5,+5.5);
	TH1D HnHitsInThisCell("HnHitsInThisCell","HnHitsInThisCell", 13, -0.5,+12.5);
	TH1D HnBestCellHits("HnBestCellHits","HnBestCellHits", 13, -0.5,+12.5);
	TH1D HCellStat("HCellStat","HCellStat",21, -0.5, 20.5);
*/	
	TH1D HFitLayers("HFitLayers","HFitLayers",21, -0.5, 20.5);
	TH1D HFitChi2("HFitChi2","HFitChi2",5000, 0., 5000.);
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
	
	
		
	
/*	
	TH1D HFitLayersAll("HFitLayersAll","HFitLayersAll",11, -0.5, 10.5);
	TH1D HFitChi2All("HFitChi2All","HFitChi2All",5000, 0., 5000.);
	TH1D HFit5Chi2All("HFit5Chi2All","HFit5Chi2All",5000, 0., 5000.);
	TH1D HFit6Chi2All("HFit6Chi2All","HFit6Chi2All",5000, 0., 5000.);
	TH1D HFit7Chi2All("HFit7Chi2All","HFit7Chi2All",5000, 0., 5000.);
	TH1D HFit8Chi2All("HFit8Chi2All","HFit8Chi2All",5000, 0., 5000.);
	
*/	
	
	
	
	TH1D HnFoundTracks("HnFoundTracks","HnFoundTracks",3, -0.5, 2.5);
	
	
	
	////////////////////////////////////////////////////////////////////
		
	// Track parameter resolutions
			
	////////////////////////////////////////////////////////////////////
	
	
	TH1D HDeltaPhi("HDeltaPhi","HDeltaPhi",100, -.005, +.005);
	TH1D HDeltaEta("HDeltaEta","HDeltaEta",100, -.005, +.005);
	TH1D HDeltaInvPt("HDeltaInvPt","HDeltaInvPt",100, -.005, +.005);
	TH1D HDeltaZ0("HDeltaZ0","HDeltaZ0",100, -1., +1.);
	TH1D HDeltaT0("HDeltaT0","HDeltaT0",100, -50., +50.);
	
	
	unsigned NchanRes = 20; // number of channels for resolutions as a function of track parameters

	
	// normalization for resolutions as a function of track parameters
	
	TH1D HPhiN("HPhiN","HPhiN", NchanRes, par.geo_gen_t_phi-par.geo_gen_t_deltaPhi, par.geo_gen_t_phi+par.geo_gen_t_deltaPhi);
	HPhiN.SetStats(false);
	TH1D HEtaN("HEtaN","HEtaN", NchanRes, par.geo_gen_t_eta-par.geo_gen_t_deltaEta, par.geo_gen_t_eta+par.geo_gen_t_deltaEta);
	HEtaN.SetStats(false);
	TH1D HInvPtN("HInvPtN","HInvPtN", NchanRes, par.geo_gen_t_invPt_min , par.geo_gen_t_invPt_max);
	HInvPtN.SetStats(false);	
		
	// bias as a function of track parameters
	
	TH1D HBiasPhiVsPhi("HBiasPhiVsPhi","HBiasPhiVsPhi", NchanRes, par.geo_gen_t_phi-par.geo_gen_t_deltaPhi, par.geo_gen_t_phi+par.geo_gen_t_deltaPhi);
	HBiasPhiVsPhi.SetStats(false);
	TH1D HBiasEtaVsPhi("HBiasEtaVsPhi","HBiasEtaVsPhi", NchanRes, par.geo_gen_t_phi-par.geo_gen_t_deltaPhi, par.geo_gen_t_phi+par.geo_gen_t_deltaPhi);
	HBiasEtaVsPhi.SetStats(false);
	TH1D HBiasInvPtVsPhi("HBiasInvPtVsPhi","HBiasInvPtVsPhi", NchanRes,  par.geo_gen_t_phi-par.geo_gen_t_deltaPhi, par.geo_gen_t_phi+par.geo_gen_t_deltaPhi);
	HBiasInvPtVsPhi.SetStats(false);
	
	TH1D HBiasPhiVsEta("HBiasPhiVsEta","HBiasPhiVsEta", NchanRes, par.geo_gen_t_eta-par.geo_gen_t_deltaEta, par.geo_gen_t_eta+par.geo_gen_t_deltaEta);
	HBiasPhiVsEta.SetStats(false);
	TH1D HBiasEtaVsEta("HBiasEtaVsEta","HBiasEtaVsEta", NchanRes, par.geo_gen_t_eta-par.geo_gen_t_deltaEta, par.geo_gen_t_eta+par.geo_gen_t_deltaEta);
	HBiasEtaVsEta.SetStats(false);
	TH1D HBiasInvPtVsEta("HBiasInvPtVsEta","HBiasInvPtVsEta", NchanRes, par.geo_gen_t_eta-par.geo_gen_t_deltaEta, par.geo_gen_t_eta+par.geo_gen_t_deltaEta);
	HBiasInvPtVsEta.SetStats(false);
		
	TH1D HBiasPhiVsInvPt("HBiasPhiVsInvPt","HBiasPhiVsInvPt", NchanRes, par.geo_gen_t_invPt_min , par.geo_gen_t_invPt_max);
	HBiasPhiVsInvPt.SetStats(false);
	TH1D HBiasEtaVsInvPt("HBiasEtaVsInvPt","HBiasEtaVsInvPt", NchanRes, par.geo_gen_t_invPt_min , par.geo_gen_t_invPt_max);
	HBiasEtaVsInvPt.SetStats(false);
	TH1D HBiasInvPtVsInvPt("HBiasInvPtVsInvPt","HBiasInvPtVsInvPt", NchanRes, par.geo_gen_t_invPt_min , par.geo_gen_t_invPt_max);
	HBiasInvPtVsInvPt.SetStats(false);	
		
	// sigma as a function of track parameters
	// temporarily holds delta squared 
	
	TH1D HSigmaPhiVsPhi("HSigmaPhiVsPhi","HSigmaPhiVsPhi", NchanRes, par.geo_gen_t_phi-par.geo_gen_t_deltaPhi, par.geo_gen_t_phi+par.geo_gen_t_deltaPhi);
	HSigmaPhiVsPhi.SetStats(false);
	TH1D HSigmaEtaVsPhi("HSigmaEtaVsPhi","HSigmaEtaVsPhi", NchanRes, par.geo_gen_t_phi-par.geo_gen_t_deltaPhi, par.geo_gen_t_phi+par.geo_gen_t_deltaPhi);
	HSigmaEtaVsPhi.SetStats(false);
	TH1D HSigmaInvPtVsPhi("HSigmaInvPtVsPhi","HSigmaInvPtVsPhi", NchanRes,  par.geo_gen_t_phi-par.geo_gen_t_deltaPhi, par.geo_gen_t_phi+par.geo_gen_t_deltaPhi);
	HSigmaInvPtVsPhi.SetStats(false);
	
	TH1D HSigmaPhiVsEta("HSigmaPhiVsEta","HSigmaPhiVsEta", NchanRes, par.geo_gen_t_eta-par.geo_gen_t_deltaEta, par.geo_gen_t_eta+par.geo_gen_t_deltaEta);
	HSigmaPhiVsEta.SetStats(false);
	TH1D HSigmaEtaVsEta("HSigmaEtaVsEta","HSigmaEtaVsEta", NchanRes, par.geo_gen_t_eta-par.geo_gen_t_deltaEta, par.geo_gen_t_eta+par.geo_gen_t_deltaEta);
	HSigmaEtaVsEta.SetStats(false);
	TH1D HSigmaInvPtVsEta("HSigmaInvPtVsEta","HSigmaInvPtVsEta", NchanRes, par.geo_gen_t_eta-par.geo_gen_t_deltaEta, par.geo_gen_t_eta+par.geo_gen_t_deltaEta);
	HSigmaInvPtVsEta.SetStats(false);
		
	TH1D HSigmaPhiVsInvPt("HSigmaPhiVsInvPt","HSigmaPhiVsInvPt", NchanRes, par.geo_gen_t_invPt_min , par.geo_gen_t_invPt_max);
	HSigmaPhiVsInvPt.SetStats(false);
	TH1D HSigmaEtaVsInvPt("HSigmaEtaVsInvPt","HSigmaEtaVsInvPt", NchanRes, par.geo_gen_t_invPt_min , par.geo_gen_t_invPt_max);
	HSigmaEtaVsInvPt.SetStats(false);
	TH1D HSigmaInvPtVsInvPt("HSigmaInvPtVsInvPt","HSigmaInvPtVsInvPt", NchanRes, par.geo_gen_t_invPt_min , par.geo_gen_t_invPt_max);
	HSigmaInvPtVsInvPt.SetStats(false);
		
	// delta squared squared for resolutions as a function of track parameters
	
	TH1D HD4PhiVsPhi("HD4PhiVsPhi","HD4PhiVsPhi", NchanRes, par.geo_gen_t_phi-par.geo_gen_t_deltaPhi, par.geo_gen_t_phi+par.geo_gen_t_deltaPhi);
	HD4PhiVsPhi.SetStats(false);
	TH1D HD4EtaVsPhi("HD4EtaVsPhi","HD4EtaVsPhi", NchanRes, par.geo_gen_t_phi-par.geo_gen_t_deltaPhi, par.geo_gen_t_phi+par.geo_gen_t_deltaPhi);
	HD4EtaVsPhi.SetStats(false);
	TH1D HD4InvPtVsPhi("HD4InvPtVsPhi","HD4InvPtVsPhi", NchanRes,  par.geo_gen_t_phi-par.geo_gen_t_deltaPhi, par.geo_gen_t_phi+par.geo_gen_t_deltaPhi);
	HD4InvPtVsPhi.SetStats(false);
	
	TH1D HD4PhiVsEta("HD4PhiVsEta","HD4PhiVsEta", NchanRes, par.geo_gen_t_eta-par.geo_gen_t_deltaEta, par.geo_gen_t_eta+par.geo_gen_t_deltaEta);
	HD4PhiVsEta.SetStats(false);
	TH1D HD4EtaVsEta("HD4EtaVsEta","HD4EtaVsEta", NchanRes, par.geo_gen_t_eta-par.geo_gen_t_deltaEta, par.geo_gen_t_eta+par.geo_gen_t_deltaEta);
	HD4EtaVsEta.SetStats(false);
	TH1D HD4InvPtVsEta("HD4InvPtVsEta","HD4InvPtVsEta", NchanRes, par.geo_gen_t_eta-par.geo_gen_t_deltaEta, par.geo_gen_t_eta+par.geo_gen_t_deltaEta);
	HD4InvPtVsEta.SetStats(false);
		
	TH1D HD4PhiVsInvPt("HD4PhiVsInvPt","HD4PhiVsInvPt", NchanRes, par.geo_gen_t_invPt_min , par.geo_gen_t_invPt_max);
	HD4PhiVsInvPt.SetStats(false);
	TH1D HD4EtaVsInvPt("HD4EtaVsInvPt","HD4EtaVsInvPt", NchanRes, par.geo_gen_t_invPt_min , par.geo_gen_t_invPt_max);
	HD4EtaVsInvPt.SetStats(false);
	TH1D HD4InvPtVsInvPt("HD4InvPtVsInvPt","HD4InvPtVsInvPt", NchanRes, par.geo_gen_t_invPt_min , par.geo_gen_t_invPt_max);
	HD4InvPtVsInvPt.SetStats(false);
	

	
	////////////////////////////////////////////////////////////////////
		
	// PT resolutions vs. abs eta and PT
			
	////////////////////////////////////////////////////////////////////
	
	
	unsigned NchanAbsRes = 20; // number of channels for PT resolutions as a function of Abs eta and PT
	
	double maxAbsEtaRes = 2.;
	//double maxAbsInvPt = 1./3.;
	double maxAbsInvPt = max(fabs(par.geo_gen_t_invPt_max), fabs(par.geo_gen_t_invPt_min));

	TH1D HAbsEtaN("HAbsEtaN","HAbsEtaN", NchanAbsRes, 0.,maxAbsEtaRes);
	HAbsEtaN.SetStats(false);
	TH1D HAbsInvPtN("HAbsInvPtN","HAbsInvPtN", NchanAbsRes, 0., maxAbsInvPt);
	HAbsInvPtN.SetStats(false);
	
	TH1D HBiasInvPtVsAbsEta("HBiasInvPtVsAbsEta","HBiasInvPtVsAbsEta", NchanAbsRes, 0., maxAbsEtaRes);
	HBiasInvPtVsAbsEta.SetStats(false);
	TH1D HBiasInvPtVsAbsInvPt("HBiasInvPtVsAbsInvPt","HBiasInvPtVsAbsInvPt", NchanAbsRes, 0., maxAbsInvPt);
	HBiasInvPtVsAbsInvPt.SetStats(false);	
	
	TH1D HSigmaInvPtVsAbsEta("HSigmaInvPtVsAbsEta","HSigmaInvPtVsAbsEta", NchanAbsRes, 0., maxAbsEtaRes);
	HSigmaInvPtVsEta.SetStats(false);
	TH1D HSigmaInvPtVsAbsInvPt("HSigmaInvPtVsAbsInvPt","HSigmaInvPtVsAbsInvPt", NchanAbsRes, 0., maxAbsInvPt);
	HSigmaInvPtVsAbsInvPt.SetStats(false);

	TH1D HD4InvPtVsAbsEta("HD4InvPtVsAbsEta","HD4InvPtVsAbsEta", NchanAbsRes,  0., maxAbsEtaRes);
	HD4InvPtVsAbsEta.SetStats(false);
	TH1D HD4InvPtVsAbsInvPt("HD4InvPtVsAbsInvPt","HD4InvPtVsAbsInvPt", NchanAbsRes, 0., maxAbsInvPt);
	HD4InvPtVsAbsInvPt.SetStats(false);
	

	
	
	
	TH1D HTrackMass("HTrackMass","HTrackMass",100,0., 1.);
	
	
	TH1D HBarrelFraction("HBarrelFraction","HBarrelFraction", 11, -0.05, 1.05);HBarrelFraction.SetStats(false);
	TH2D HNDvsNB("NDvsNB","NDvsNB",  13, -0.5,+12.5, 13, -0.5, +12.5); HNDvsNB.SetStats(false);
	
	TH1D HlayerInd("layerInd","layerInd",  51, -0.5,+50.5); HlayerInd.SetStats(false);
	
	
	
	TH1D HitX("HitX","HitX", 4000, -2000.,+2000.); HitX.SetStats(false);
	TH1D HitY("HitY","HitY", 4000, -2000.,+2000.); HitY.SetStats(false);
	TH1D HitZ("HitZ","HitZ", 5000, -2500.,+2500.); HitZ.SetStats(false);
	
	// Track hits
	TH2D HitTXYbarrel("HitTXYbarrel","HitTXYbarrel", 4000, -2000.,+2000., 4000, -2000, +2000); HitTXYbarrel.SetStats(false);
	TH2D HitTXYdisc("HitTXYdisc","HitTXYdisc", 4000, -2000.,+2000., 4000, -2000, +2000); HitTXYdisc.SetStats(false);
	// Bib hits
	TH2D HitBXYbarrel("HitBXYbarrel","HitBXYbarrel", 4000, -2000.,+2000., 4000, -2000, +2000); HitBXYbarrel.SetStats(false);
	TH2D HitBXYdisc("HitBXYdisc","HitBXYdisc", 4000, -2000.,+2000., 4000, -2000, +2000); HitBXYdisc.SetStats(false);
	
	
	
	TH2D HitTRZ("HitTRZ","HitTRZ",  5000, -2500.,+2500., 1800, 0, 1800); HitTRZ.SetStats(false);
	TH2D HitBRZ("HitBRZ","HitBRZ",  5000, -2500.,+2500., 1800, 0, 1800); HitBRZ.SetStats(false);
	
	
	TH1D HitPhi("HitPhi","HitPhi", 1000, -Pi,+Pi);HitPhi.SetStats(false);
	TH2I HitZPhi("HitZPhi","HitZPhi",  10000, -2500.,+2500., 10000, -Pi, +Pi); HitZPhi.SetStats(false);
	TH2I HitRPhi("HitRPhi","HitRPhi",  10000, 0.,+2500., 10000, -1., +1.); HitRPhi.SetStats(false);
	
	
	TH1D HTrackZ0("HTrackZ0","HTrackZ0", 600, -300,+300); HTrackZ0.SetStats(true);
	TH1D HTrackT0("HTrackT0","HTrackT0", 600, -300,+300); HTrackT0.SetStats(true);
	
	TH1D HTrackP("HTrackP","HTrackP", 1000, 0.,10.); HTrackP.SetStats(true);
	
	
	
	unsigned NchanEff = 20; // number of channels for efficiency as a function of track parameters
	
	
	TH1D HTrackEta("HTrackEta","HTrackEta", NchanEff, par.geo_gen_t_eta-par.geo_gen_t_deltaEta, par.geo_gen_t_eta+par.geo_gen_t_deltaEta);
	TH1D HTrackPhi("HTrackPhi","HTrackPhi", NchanEff, par.geo_gen_t_phi-par.geo_gen_t_deltaPhi, par.geo_gen_t_phi+par.geo_gen_t_deltaPhi);
	TH1D HTrackInvPt("HTrackInvPt","HTrackInvPt", NchanEff, par.geo_gen_t_invPt_min , par.geo_gen_t_invPt_max);

	TH1D HTrackNhitsVsEta("HTrackNhitsVsEta","HTrackNhitsVsEta", NchanEff, par.geo_gen_t_eta-par.geo_gen_t_deltaEta, par.geo_gen_t_eta+par.geo_gen_t_deltaEta);
	TH1D HTrackNhitsVsInvPt("HTrackNhitsVsInvPt","HTrackNhitsVsInvPt", NchanEff, par.geo_gen_t_invPt_min , par.geo_gen_t_invPt_max);
	TH1D HTrackNhits2VsEta("HTrackNhits2VsEta","HTrackNhits2VsEta", NchanEff, par.geo_gen_t_eta-par.geo_gen_t_deltaEta, par.geo_gen_t_eta+par.geo_gen_t_deltaEta);
	TH1D HTrackNhits2VsInvPt("HTrackNhits2VsInvPt","HTrackNhits2VsInvPt", NchanEff, par.geo_gen_t_invPt_min , par.geo_gen_t_invPt_max);
	
	
	TH1D HTrackEtaEff("HTrackEtaEff","HTrackEtaEff", NchanEff, par.geo_gen_t_eta-par.geo_gen_t_deltaEta, par.geo_gen_t_eta+par.geo_gen_t_deltaEta);
	TH1D HTrackPhiEff("HTrackPhiEff","HTrackPhiEff", NchanEff, par.geo_gen_t_phi-par.geo_gen_t_deltaPhi, par.geo_gen_t_phi+par.geo_gen_t_deltaPhi);
	TH1D HTrackInvPtEff("HTrackInvPtEff","HTrackInvPtEff", NchanEff, par.geo_gen_t_invPt_min , par.geo_gen_t_invPt_max);
	
	
		
	unsigned NchanAbsEff = 20; // number of channels for efficiency as a function of Abs eta and PT
	//double maxAbsEtaEff = 2.5; // max abs eta for efficiency plots
	double maxAbsEtaEff = fabs(par.geo_gen_t_eta + par.geo_gen_t_deltaEta);// max abs eta for efficiency plots
	
		
	TH1D HTrackAbsEta("HTrackAbsEta","HTrackAbsEta", NchanAbsEff,0.,maxAbsEtaEff); 
	TH1D HTrackAbsInvPt("HTrackAbsInvPt","HTrackAbsInvPt", NchanAbsEff, 0., maxAbsInvPt);

	TH1D HTrackAbsEtaEff("HTrackAbsEtaEff","HTrackAbsEtaEff", NchanAbsEff ,0.,maxAbsEtaEff); 
	TH1D HTrackAbsInvPtEff("HTrackAbsInvPtEff","HTrackAbsInvPtEff", NchanAbsEff, 0., maxAbsInvPt);

	// fraction of BIB hits vs eta and PT
	
	TH1D HHitAbsEta("HHitAbsEta","HHitAbsEta", NchanAbsEff,0.,maxAbsEtaEff); 
	TH1D HHitAbsInvPt("HHitAbsInvPt","HHitAbsInvPt", NchanAbsEff, 0., maxAbsInvPt);

	TH1D HHitAbsEtaBIB("HHitAbsEtaBIB","HHitAbsEtaBIB", NchanAbsEff ,0.,maxAbsEtaEff); 
	TH1D HHitAbsInvPtBIB("HHitAbsInvPtBIB","HHitAbsInvPtBIB", NchanAbsEff, 0., maxAbsInvPt);



	
	TH2I HTrackInvPtvsPhi("HTrackInvPtvsPhi","HTrackInvPtvsPhi",1000, -Pi,+Pi,1000, -1/2.,+1/2.);HTrackInvPtvsPhi.SetStats(true);
	
	TH1D HTrackPz("HTrackPz","HTrackPz", 10000, -100., +100.); HTrackPz.SetStats(true);
	TH1D HTrackInvPz("HTrackInvPz","HTrackInvPz", 10000, -5., +5.); HTrackInvPz.SetStats(true);
	TH1D HTrackPhi0("HTrackPhi0","HTrackPhi0", 1000, -Pi,+Pi); HTrackPhi0.SetStats(true);
	TH2I HTrackInvPzvsPhi0("HTrackInvPzvsPhi0","HTrackInvPzvsPhi0",1000, -Pi,+Pi,1000, -1/2.,+1/2.);HTrackInvPzvsPhi0.SetStats(true);
	
	TH1D HTrackLayersHit("HLayersHit","HLayersHit", 13, -0.5,12.5);
	TH2I HTrackLayersHitVsTrackLayers("HLayersHitVsLayers","HLayersHitVsLayers", 13, -0.5,12.5, 13, -0.5,12.5);  
			HTrackLayersHitVsTrackLayers.SetStats(false);
	TH2I HTrackLayersHitVsTrackHits("HLayersHitVsTrackHits","HLayersHitVsTrackHits", 13, -0.5,12.5, 13, -0.5,12.5);  
			HTrackLayersHitVsTrackLayers.SetStats(false);
	TH2I HTCellLayersHitVsMinLayers("HTCellLayersHitVsMinLayers","HTCellLayersHitVsMinLayers", 13, -0.5,12.5, 13, -0.5,12.5);  
			HTCellLayersHitVsMinLayers.SetStats(false);
	TH1D HTMissingLayers("HTMissingLayers","HTMissingLayers",18, -3.5,14.5);
	
	
	TH2I H2nHitVSnMin("H2nHitVSnMin","H2nHitVSnMin",  21, -0.5,+20.5, 21, -0.5,+20.5);
	
	TH1D HDistanceToBestCell("HDistanceToBestCell","HDistanceToBestCell",11, -0.5,10.5);


/*
	
	// Define vector of 1-Dim histograms for barrel hits
	
	const static unsigned nBarrels = dg.nBarrels;
	
	std::vector<TH1D*> HitTBX; // Track hits	
	for(unsigned iB = 0; iB != nBarrels; ++iB){
		stringstream ss;
		ss << "HitTBX" << iB;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, 2* dg.B[iB].r * Pi, -dg.B[iB].r * Pi, dg.B[iB].r * Pi); 
		HitTBX.push_back(h);	
	}
	
	std::vector<TH1D*> HitBBX; // BIB hits	
	for(unsigned iB = 0; iB != nBarrels; ++iB){
		stringstream ss;
		ss << "HitBBX" << iB;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, 2* dg.B[iB].r * Pi, -dg.B[iB].r * Pi, dg.B[iB].r * Pi); 
		HitBBX.push_back(h);	
	}
	
	
	
	std::vector<TH1D*> HitTBZ; // Track Hits	
	for(unsigned iB = 0; iB != nBarrels; ++iB){
		stringstream ss;
		ss << "HitTBZ" << iB;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, -dg.B[iB].zMin + dg.B[iB].zMax, dg.B[iB].zMin, dg.B[iB].zMax); 
		HitTBZ.push_back(h);	
	}
	
	std::vector<TH1D*> HitBBZ; // Bib Hits	
	for(unsigned iB = 0; iB != nBarrels; ++iB){
		stringstream ss;
		ss << "HitBBZ" << iB;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, -dg.B[iB].zMin + dg.B[iB].zMax, dg.B[iB].zMin, dg.B[iB].zMax); 
		HitBBZ.push_back(h);	
	}
	
	
	
	std::vector<TH1D*> HitTBT; // Track hits	
	for(unsigned iB = 0; iB != nBarrels; ++iB){
		stringstream ss;
		ss << "HitTBT" << iB;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, 2000, 0., 2000.); 
		HitTBT.push_back(h);	
	}
	
	std::vector<TH1D*> HitBBT; // Bib hits	
	for(unsigned iB = 0; iB != nBarrels; ++iB){
		stringstream ss;
		ss << "HitBBT" << iB;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, 2000, 0., 2000.); 
		HitBBT.push_back(h);	
	}
	
	std::vector<TH1D*> HitBTX;	
	for(unsigned iB = 0; iB != nBarrels; ++iB){
		stringstream ss;
		ss << "HitBTX" << iB;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, 4000, -2000., +2000.); 
		HitBTX.push_back(h);	
	}
	
	
	
	std::vector<TH1D*> HitBXCorr;
	for(unsigned iB = 0; iB != nBarrels; ++iB){
		stringstream ss;
		ss << "HitBXCorr" << iB;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, 2* dg.B[iB].r * Pi, -dg.B[iB].r * Pi, dg.B[iB].r * Pi); 
		HitBXCorr.push_back(h);	
	}
	
	
	// Define vector of 1-Dim histograms for disc hits
	
	const static unsigned nDiscs = dg.nDiscs;
		
	std::vector<TH1D*> HitTDPhi; // track Hits
	
	for(unsigned iD = 0; iD != nDiscs; ++iD){
		stringstream ss;
		ss << "HitTDPhi" << iD;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, 1000, -1., +1.);
		HitTDPhi.push_back(h);	
	}
	
	std::vector<TH1D*> HitBDPhi; // Bib hits
	
	for(unsigned iD = 0; iD != nDiscs; ++iD){
		stringstream ss;
		ss << "HitBDPhi" << iD;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, 1000, -1., +1.);
		HitBDPhi.push_back(h);	
	}
	
	std::vector<TH1D*> HitTDR; // Track hits
	
	for(unsigned iD = 0; iD != nDiscs; ++iD){
		stringstream ss;
		ss << "HitTDR" << iD;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, -dg.D[iD].rMin + dg.D[iD].rMax, dg.D[iD].rMin, dg.D[iD].rMax);
		HitTDR.push_back(h);	
	}
	
	std::vector<TH1D*> HitBDR; // Bib hits
	
	for(unsigned iD = 0; iD != nDiscs; ++iD){
		stringstream ss;
		ss << "HitBDR" << iD;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, -dg.D[iD].rMin + dg.D[iD].rMax, dg.D[iD].rMin, dg.D[iD].rMax);
		HitBDR.push_back(h);	
	}
	
		
	std::vector<TH1D*> HitTDT; // Track hits	
	for(unsigned iD = 0; iD != nDiscs; ++iD){
		stringstream ss;
		ss << "HitTDT" << iD;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, 4000, 0., 4000.); 
		HitTDT.push_back(h);	
	}
			
	std::vector<TH1D*> HitBDT; // Bib hits	
	for(unsigned iD = 0; iD != nDiscs; ++iD){
		stringstream ss;
		ss << "HitBDT" << iD;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, 4000, 0., 4000.); 
		HitBDT.push_back(h);	
	}
	
		
	
	std::vector<TH1D*> HitDTX;	
	for(unsigned iD = 0; iD != nDiscs; ++iD){
		stringstream ss;
		ss << "HitDTX" << iD;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, 4000, -2000., 2000.); 
		HitDTX.push_back(h);	
	}
	
	
			
	std::vector<TH1D*> HitDPhiCorr;
	
	for(unsigned iD = 0; iD != nDiscs; ++iD){
		stringstream ss;
		ss << "HitDPhiCorr" << iD;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, 1000, -1., +1.);
		HitDPhiCorr.push_back(h);	
	}
	
	
	
	// Define vector of 2-Dim histograms for barrel hits
	
	
	std::vector<TH2I*> HitBXZ;
	
	for(unsigned iB = 0; iB != nBarrels; ++iB){
		stringstream ss;
		ss << "HitBXZ" << iB;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH2I* h = new TH2I(title,title, 2* dg.B[iB].r * Pi, -dg.B[iB].r * Pi, dg.B[iB].r * Pi,1000,dg.B[iB].zMin,dg.B[iB].zMax);
		h->SetStats(false);
		HitBXZ.push_back(h);	
	}
	
	// Define vector of 2-Dim histograms for disc hits
	
	std::vector<TH2I*> HitDPhiR;
	
	for(unsigned iD = 0; iD != nDiscs; ++iD){
		stringstream ss;
		ss << "HitDPhiR" << iD;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH2I* h = new TH2I(title,title,1000, -1., +1., 1000,dg.D[iD].rMin,dg.D[iD].rMax);
		h->SetStats(false);
		HitDPhiR.push_back(h);	
	}
	
*/
	

	
	// MAIN LOOP ON EVENT GENERATION  ///////////////////////////////////////////////
	
	
	struct FitTrack {
	
        		double chi2;
        		double phi;
        		double eta;
        		double invPt; 
        		double z0;
        		double t0;
        		double mass; 
        		unsigned nLayers;
    };
  	
	
	
	double candidateRate = 0.;
	unsigned nEventsWithCandidates = 0;
	unsigned nEventsWithTracks = 0;
	
	// instantiate track reader from file
	
	TrackReader* newRecoTrack = new TrackReader(par.reco_inputTrackFileName);
	
	
	int argc; 
	char **argv;

//	TApplication theApp("app",&argc, argv);
	
	TGraph2D* g2;
  	TCanvas* c;
 
	
	cout << "BEGIN EVENT GENERATION " << endl;	

    for(unsigned iEv = 0; iEv != nEvents; ++iEv){
    
    	
	cout << "********************************************************************" << endl;	
		
   			// declaring argument of time()
    		time_t my_time = time(NULL);
  
   			// ctime() used to give the present time
    		printf("%s", ctime(&my_time));
    				
    		cout << iEv << "/" << nEvents << " processed events" << endl;
    		
    		//if(nEvents >= 1e6)
    		//if(iEv && ((iEv % (int)100) == 0) ) cout << iEv << "/" << nEvents << " processed events" << endl;
    		
    	
    		// create one event
    		
    		
			Event ev(newRecoTrack, nTracks, bibRead, backGnd);// (track reader, nTracks, BIB reader, N bkg hits)
			if((int)ev.trackList.size() < (int)nTracks) break; // end of input file
			
			if(verbose) ev.print(cout,1);	
									 
			if(verbose) cout << "event generated" << endl;
    		vector <FitTrack> foundTracks; // tracks found in this event (now empty)
    		
	
			HTA->reset(); // clear all hits from HT Array
			
			if(verbose) cout << "HTA reset" << endl;
	
		// LOOP ON TRACKS 
		
		for(int iT = 0; iT != (int)ev.trackList.size(); ++iT){
		
			Track thisTrack = ev.trackList[iT];
			
			
			int i,j,k;
			trackToCell(i, j, k, ev.trackList[iT]);
			if(verbose) cout << i << " " << j << " " << k << endl;
			
			if(verbose) thisTrack.print(cout,1); // print track 
			
			// Fill histograms of main track parameters	
			
			HTrackZ0.Fill(thisTrack.z0);
			HTrackT0.Fill(thisTrack.t0);
			HTrackEta.Fill(thisTrack.eta);
			HTrackAbsEta.Fill(fabs(thisTrack.eta));
			HTrackPhi.Fill(thisTrack.phi);
			HTrackInvPt.Fill(thisTrack.invPt);
			HTrackAbsInvPt.Fill(fabs(thisTrack.invPt));
			HTrackInvPtvsPhi.Fill(thisTrack.phi,thisTrack.invPt);
			
			HTrackPz.Fill(thisTrack.Pz);
			HTrackInvPz.Fill(1./thisTrack.Pz);
			HTrackPhi0.Fill(thisTrack.phi0);
			HTrackInvPzvsPhi0.Fill(thisTrack.phi0,1./thisTrack.Pz);
			
			// average number of hits vs. eta and pt
			
			unsigned numHits = thisTrack.hitList.size();
				HTrackNhitsVsEta.Fill(thisTrack.eta,(double)numHits);
				HTrackNhits2VsEta.Fill(thisTrack.eta,(double)(numHits*numHits));
				HTrackNhitsVsInvPt.Fill(thisTrack.invPt,(double)numHits);
				HTrackNhits2VsInvPt.Fill(thisTrack.invPt,(double)(numHits*numHits));					
			
		}
	
		
		// LOOP ON HITS
	
		unsigned nHits = ev.hitList.size();
		
		if(verbose) cout << "nHits " << nHits << endl;
		
		for(unsigned iH = 0; iH != nHits; ++iH) {
			
			Hit thisHit = ev.hitList[iH];
			if(!thisHit.isSeed()) continue;// ignore all non-seed hits
			
		
	/*		if(par.reco_fillHitHistograms){ // fill all hit histograms
			
				double X, Y, Z;
				thisHit.XYZ(dg, X, Y, Z);
				double R = sqrt(X*X + Y*Y);
				double Phi = atan2(Y,X);
				
				HlayerInd.Fill(thisHit.layerInd);
			
			
				HitX.Fill(X);
				HitY.Fill(Y);
				HitZ.Fill(Z);
				HitPhi.Fill(Phi);
				HitZPhi.Fill(Z,Phi);
				HitRPhi.Fill(R,Phi);
			
				if(thisHit.hitType == 'B') { 
				
					HitBXZ[thisHit.iLayer]->Fill(thisHit.x1,Z);
					if(thisHit.trackInd == 0) {
						HitBRZ.Fill(Z,R);
						HitBXYbarrel.Fill(X,Y);
						HitBBX[thisHit.iLayer]->Fill(thisHit.x1);
						HitBBZ[thisHit.iLayer]->Fill(thisHit.x2);
						HitBBT[thisHit.iLayer]->Fill(thisHit.t);
					}
					else {
						HitTRZ.Fill(Z,R);
						HitTXYbarrel.Fill(X,Y);
						HitTBX[thisHit.iLayer]->Fill(thisHit.x1);
						HitTBZ[thisHit.iLayer]->Fill(thisHit.x2);
						HitTBT[thisHit.iLayer]->Fill(thisHit.t);
					}
				}
				if(thisHit.hitType == 'D') {
					HitDPhiR[thisHit.iLayer]->Fill(Phi,R);	
					if(thisHit.trackInd == 0) {				
						HitBRZ.Fill(Z,R);		
						HitBDPhi[thisHit.iLayer]->Fill(Phi);	
						HitBXYdisc.Fill(X,Y);
						HitBDR[thisHit.iLayer]->Fill(R);
						HitBDT[thisHit.iLayer]->Fill(thisHit.t);
						HitBDPhi[thisHit.iLayer]->Fill(Phi);		
					}
					else {				
						HitTRZ.Fill(Z,R);		
						HitTDPhi[thisHit.iLayer]->Fill(Phi);		
						HitTXYdisc.Fill(X,Y);
						HitTDR[thisHit.iLayer]->Fill(R);
						HitTDT[thisHit.iLayer]->Fill(thisHit.t);
						HitTDPhi[thisHit.iLayer]->Fill(Phi);					
					}
				}
			} // end fill all hit histograms
			
	*/		//////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////			
			//////////////////////////////////////////////////////////////////////////////
			/// FILL HTA ARRAY FOR PATTERN REGOGNITION ///////////////////////////////////
		
			
			if(thisHit.isSeed())HTA->fill(thisHit, par.reco_fillMode);// use only seed hits 
		
			Special = false;
			
				
		}// end loop on hits
		
		
		
		if(Debug) cout << endl;
		
		unsigned phi_b, eta_b, pt_b;// best cell
		
		unsigned nLayersHit = HTA->getBestCell(phi_b, eta_b, pt_b);// num of layers hit in best cell
		
		if(true) cout << endl << "bestCell: [" << phi_b << "," << eta_b << "," << pt_b << "] nLayersHit: " << nLayersHit << endl;	
		//HTA->ArrElem[phi_b][eta_b][pt_b].printHits(cout);
		
		
		if(nLayersHit == 0) continue; // skip to next event		
		//HnBestCellHits.Fill(nLayersHit);
		
		 	bool foundTrack = false;
		 	
		 	//////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////
			/////// CALL fitCandidate FOR BEST CELL
							
			if(par.reco_TrackFit){
										
				double chi2, phi, eta, invPt, z0, t0, beta;
				vector <long int> goodFitHitList;
				unsigned nHitsFit;
					
				//cout << "FitCandidate in"	<< 	endl;
						
				int retCodeFit = HTA->ArrElem[phi_b][eta_b][pt_b].
						fitCandidate(&(HTA->allHits),chi2, phi, eta, invPt, z0, t0, beta, nHitsFit, goodFitHitList);
				
				//cout << "FitCandidate out"	<< 	endl;
						
				unsigned dof = 2*nHitsFit - 5; // Degrees of freedom of the fit
				
				cout << endl << "retCodeFit: " << retCodeFit << " chi2 = "<< chi2 << " d.o.f. = " << dof << endl;
				HFitChi2.Fill(chi2);
				
				Track thisTrack = ev.trackList[0];
				cout << endl << "Generated track: phi = " << thisTrack.phi << " eta = " << thisTrack.eta << " invPt = " << thisTrack.invPt 
					<< " z0 = " << thisTrack.z0 << " t0 = " << thisTrack.t0 << " beta = " << thisTrack.beta << endl;
										
				cout << "Fit result:      phi = " << phi << " eta = " << eta << " invPt = " << invPt 
					<< " z0 = " << z0 << " t0 = " << t0 << " beta = " << beta << endl;
		
		
				if(retCodeFit != 0){
		
					cout << "FIT FAILED with retCodeFit = " << retCodeFit << endl;			
					continue;		
				}
		
					
				if(par.reco_plotTrack3D){	
					
					/////////////////////////////////////////////////////////////////	
					///// display track and pause
				
					static const unsigned MAXLAYERS = 100;
				
					double xp[MAXLAYERS], yp[MAXLAYERS],zp[MAXLAYERS];			

					for(int ih = 0; ih != nHitsFit; ++ih){// begin loop on hits				  
					  Hit* thisHit = &ev.hitList[goodFitHitList[ih]];				  
					  xp[ih] = thisHit->x;
					  yp[ih] = thisHit->y;
					  zp[ih] = thisHit->z;
				  
					}// end loop on hits
				

					c = new TCanvas("c","Track plot",0,0,800,1250);

					g2 = new TGraph2D(nHitsFit, xp, yp, zp);

					g2->SetMarkerStyle(20);
					g2->SetMarkerSize(1);

					TH3I *Track_Plot = new TH3I("Track_Plot","Track_Plot",1,-1600.,+1600.,1,-1600.,1600.,1,-2500.,+2500.);
					Track_Plot->SetStats(false);

					Track_Plot->SetTitle("; X [mm];Y [mm];Z [mm]");

					Track_Plot->Draw();

					g2->Draw("PSAME");

					gPad->Update();
					gPad->WaitPrimitive();

					delete c;
					delete g2;
					delete Track_Plot;

				  
				} // end if(par.reco_plotTrack3D)
				
				  
				////////////////////////////////////////////////////////////////////  
				// FILL HISTOGRAMS
				  
				  // TRACK PARAMETER RESOLUTIONS
				  		  
					HDeltaPhi.Fill(phi-thisTrack.phi );
					HDeltaEta.Fill(eta-thisTrack.eta);
					HDeltaInvPt.Fill(invPt-thisTrack.invPt);
					HDeltaZ0.Fill(z0-thisTrack.z0);
					HDeltaT0.Fill(t0-thisTrack.t0);
	


			}// end if(par.reco_TrackFit)							
											
	
	} // end loop on events ///////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////
	
	
 	
  
	
	cout << "***************************************************" << endl;
	cout << "EVENT GENERATION COMPLETE" << endl;
	
	
			
	/////  Efficiency ///////////////////////////////////////
	
	cout << nEventsWithTracks << " events with tracks out of " << nEvents << endl;
	cout << "Efficiency = " << nEventsWithTracks/(double)nEvents << endl;

			
	///// Events ///////////////////////////////////
	
	cout << "Number of Events with candidates = " << nEventsWithCandidates << endl;
	//cout << "Trigger probability = " << nEventsWithCandidates/(double)nEvents << endl;
			
	///// Rate (candidates/collision) /////////////////////////////////////////////
			
	cout << "Total number of candidates = " << candidateRate << endl;
	//cout << "Number of candidates per Collision = " << candidateRate/nEvents*rateScale << endl;
	
	
	

		
	cout << "Writing histogram file " << histFileName << "..." << endl;
	histFile->Write(); // write histogram file
	

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

	// theApp.Run();
	
	return 0;		
	
} // end main



