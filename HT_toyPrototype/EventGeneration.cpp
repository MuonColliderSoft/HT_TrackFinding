// simulation of track reconstruction. Aug 29 2022
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

#include "Statistics.cpp"

#include "makeEffGraphFromHists.cpp"


using namespace std;

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
//std::normal_distribution<double> gauss(0.,1.);


#include "Gauss.h"


#include "GlobalConstants.cpp"

#include "Parameters.h"
Parameters par; // class that holds all configurable parameters					   

#include "DetectorGeometry.cpp"
DetectorGeometry dg;

#include "TrackGeometry.cpp"
TrackGeometry tg("gen"); // track parameters for event generation

#include "Hit.h"
#include "Track.cpp"
#include "Hit.cpp"
#include "BibFileReader.cpp"
#include "Event.cpp"
#include "HTArray.cpp"


bool Debug = false;
bool verbose = par.gen_verbose;

int TrainingPhase = 0; // 0 = pattern recognition 
                       // 1 = from scratch
                       // 2 = second step after diagonalization
                       // 3 = verification and statistics
				   
// The following global parameters are initialized in the main function

unsigned nEvents; // Number of events to be generated
unsigned nTracks; // Number of tracks per event
double fracBib; // fraction of Bib to be simulated
int fillMode; // optimization mode for HTA fill
bool PlotTracks; // create a file with data for 3-D plots of candidates		


				   

int main(){

	
	dg.print(cout);

	// parameter initialization
    nEvents = par.gen_nEvents; // Number of events to be generated
	nTracks = par.gen_nTracks; // Number of tracks per event
	fracBib = par.gen_fracBib; // fraction of Bib to be simulated
	fillMode = par.gen_fillMode; // optimization mode for HTA fill
	PlotTracks = par.gen_PlotTracks; // create a file with data for 3-D plots of candidates
	
	// seeding the random generators 
  	long int randomSeed = par.gen_randomSeed;
  	if(randomSeed)generator.seed(randomSeed);
  	long int randomSeed_trk = par.gen_randomSeed_trk;
  	if(randomSeed_trk)generator_trk.seed(randomSeed_trk);  
  	
  	string histFileName = par.gen_histFileName;
	TFile* histFile = new TFile(histFileName.c_str(),"RECREATE");  // histogram file
  
	// Instantiate Hough Transform Array
	
	cout << "Instantiating Hough Transform Array..." << endl;
	
	static HTArray HTA;	
	HTA.initHists();
	HTA.print(cout);
			
	// Open data file
	
	string dataFileName = par.gen_dataFileName;
	cout << "Opening Data File: " << dataFileName << " ..." << endl;
	ifstream infile;
	infile.open(dataFileName);
	
	if(infile){ // data file found
		int retcode = HTA.read(infile);
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
		
	
		ofstream outPlotFile; // file for to plot candidates
		
		if(PlotTracks) {	
	
	// Open plot tracks data file
	
		string plotDataFileName = par.gen_plotDataFileName;
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
	
	if(fracBib > 0){
	
		// opening BIB background file
	
		cout << "reading BIB... ";
	
		string bibFileName = par.gen_bibFileName;	
		int code = bibRead.readFile(bibFileName);
		cout << "return code = " << code << endl;
		if(code == -1 && fracBib > 0.) return 0;
	
		unsigned totBIB = bibRead.size();
		cout << totBIB << " BIB events read from file" << endl;	
	
		// calculating the number of Bib hits to add to every event
		// given the slice in phi that we want to cover and what
		// fraction of background we want to simulate (fracBib)
	
		double phia = par.gen_phia; // low limit in phi
		double phib = par.gen_phib; // high limit in phi
		bibRead.setPhiLimits(phia,phib); 
		cout << "BIB phi limits: " << phia << " " << phib << endl;
		double BIBscale = 2.*Pi/(phib - phia);
		cout << "BIB scale factor: " << BIBscale << endl;	
		nBibHits = totBIB/BIBscale*fracBib;
			if(fracBib < 0.) nBibHits = 0;
		cout << nBibHits << " BIB scaled hits per collision" << endl;
	
	}
					
	std::poisson_distribution<unsigned> poissDist(nBibHits);
	
	
	double rateScale = 2.*Pi/(HTA.phiStep*HTA.NphiBins);
	cout << "rate scale factor: " << rateScale << endl;	
	
	
			
	TH1D HnCandidates("HnCandidates","HnCandidates", 6, -0.5,+5.5);
	TH1D HnHitsInThisCell("HnHitsInThisCell","HnHitsInThisCell", 13, -0.5,+12.5);
	TH1D HnBestCellHits("HnBestCellHits","HnBestCellHits", 13, -0.5,+12.5);
	TH1D HCellStat("HCellStat","HCellStat",21, -0.5, 20.5);
	
	TH1D HFitLayers("HFitLayers","HFitLayers",11, -0.5, 10.5);
	TH1D HFitChi2("HFitChi2","HFitChi2",5000, 0., 5000.);
	TH1D HFit5Chi2("HFit5Chi2","HFit5Chi2",5000, 0., 5000.);
	TH1D HFit6Chi2("HFit6Chi2","HFit6Chi2",5000, 0., 5000.);
	TH1D HFit7Chi2("HFit7Chi2","HFit7Chi2",5000, 0., 5000.);
	TH1D HFit8Chi2("HFit8Chi2","HFit8Chi2",5000, 0., 5000.);
	
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
	

	/*
		// cell centers
	
		cout << endl;
		for(int i = 0; i != HTA_NphiBins; ++i) cout << "i: " << i << " phi: " << HTA.getCell_Phi(i) << endl;
		cout << endl;
		for(int j = 0; j != HTA_NetaBins; ++j) cout << "j: " << j << " eta: " << HTA.getCell_Eta(j) << endl;
		cout << endl;
		for(int k = 0; k != HTA_NinvptBins; ++k)cout << "k: " << k << " invPt: " << HTA.getCell_InvPt(k) << endl;
		cout << endl;										
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
    		
			Event ev(dg, nTracks);// (geometry, nTracks)
			
			if(nBibHits){
				
				// fluctuate BIB hits			
				unsigned nxBibHits = poissDist(generator);
							
				//add BIB hits		
				ev.addBibHits(bibRead, nxBibHits);
			}
			
			 
			 
    		vector <FitTrack> foundTracks; // tracks found in this event (now empty)
	
			HTA.reset(); // clear all hits from HT Array
	
		// FORMER LOOP ON TRACKS (now only one track)
		
		if(!ev.trackList.empty()){
		
			//unsigned nTracks = 1;
			unsigned iT = 0;
			
			Track thisTrack = ev.trackList[iT];
			
			thisTrack.print(cout,1); // print track with all its hits
			
			// Fill histograms of main track parameters	
			
			HTrackZ0.Fill(thisTrack.z0);
			HTrackT0.Fill(thisTrack.t0);
			HTrackEta.Fill(thisTrack.eta);
			HTrackPhi.Fill(thisTrack.phi);
			HTrackInvPt.Fill(thisTrack.invPt);
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
		
		if(Debug) cout << "nHits " << nHits << endl;
		for(unsigned iH = 0; iH != nHits; ++iH) {
			
			Hit thisHit = ev.hitList[iH];
			if(!thisHit.isSeed()) continue;// ignore all vertex hits
			
			if(par.gen_fillHitHistograms){ // fill all hit histograms
			
				if(Debug) cout << " processing layer " << thisHit.iLayer << endl;
			
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
			
			
			/// FILL HTA ARRAY FOR PATTERN REGOGNITION ///////////////////////////////////
			
			if(thisHit.isSeed())HTA.fill(thisHit,fillMode);// use only non-vertex hits (seed hits)
			Special = false;
			
				
		}// end loop on hits
		
		
		
		if(Debug) cout << endl;
		
		unsigned phi_b, eta_b, pt_b;// best cell
		
		unsigned nLayersHit = HTA.getBestCell(phi_b, eta_b, pt_b);// num of layers hit in best cell
		if(nLayersHit == 0) continue; // skip to next event

		unsigned nLayers = HTA.ArrElem[phi_b][eta_b][pt_b].layerIndHitStat.size();// total num of layers in best cell
		unsigned minLayers = HTA.ArrElem[phi_b][eta_b][pt_b].minLayers;// min possible layers hit in best cell
		
		HnBestCellHits.Fill(nLayersHit);
		
		
		
		HTrackLayersHit.Fill(nLayersHit);// num of layers hit in best cell
		HTrackLayersHitVsTrackLayers.Fill(nLayers, nLayersHit);
		HTrackLayersHitVsTrackHits.Fill(nHits, nLayersHit);/////////////
		HTCellLayersHitVsMinLayers.Fill(minLayers, nLayersHit);
		int diff = minLayers - nLayersHit;
		HTMissingLayers.Fill(diff);
		
	
		
		///////////////////////////////////////////////////////////////
		/////////////// SUMMARY /////////////////////////////////////////
		
		//if(nHits != nLayersHit){
		
		unsigned nHitsInThisCell = 0;
		unsigned nMinLayers = 0;
		
				
			//cout << endl;
			if(verbose) cout << "bestCell " << nLayersHit << " layers hit in [" << phi_b << "][" << eta_b << "][" << pt_b << "]" << endl;	
			if(!ev.trackList.empty()){// event contains one original track - report			
				//cout << " Event " << iEv << " nHits " << nHits << " foundHits " << nLayersHit << endl;
				int it, jt, kt;
				ev.trackList[0];
				int retcode = HTA.getCell(ev.trackList[0],it,jt,kt);			
				unsigned nHitsInThisCell = 	HTA.ArrElem[it][jt][kt].nHitLayers;
				HnHitsInThisCell.Fill(nHitsInThisCell);	
				unsigned nMinLayers = HTA.ArrElem[it][jt][kt].minLayers;
					
				//cout << " track cell " << " phi " << it << " eta " << jt << " pt " << kt;
				if(verbose) cout << "Track cell " << " [" << it << "][" << jt << "][" << kt << "]";			
				if(verbose) cout << " nHitsInThisCell " << nHitsInThisCell << endl;
						
						
		
			} // end report for fate of original track	
			
			
			// Deal with candidates
							
			unsigned nCan = HTA.getCellCandidates();
			HnCandidates.Fill(nCan);			
			if(nCan) ++nEventsWithCandidates;
			
			if(verbose) cout << "Event " << iEv << ": "<< nCan << " Candidates" << endl;
					
			if(verbose) if(par.gen_printCandidates) HTA.printCellCandidateList(cout);
			
			if(PlotTracks) {
				// write file with track candidates
				outPlotFile << "Event " << iEv << " ";
				HTA.writeCellCandidateList(outPlotFile, ev.hitList);
			}
			
			cout << "Loop on candidates" << endl;
				
			bool foundTrack = false;
			
			// Loop on candidates 
			
			for(unsigned iCan = 0; iCan != nCan; ++iCan){
				HTArray::Pars p = HTA.cellCandidateList[iCan];
				unsigned minLayers = HTA.ArrElem[p.iPhi][p.iEta][p.iInvpt].minLayers;
				unsigned nHitLayers = HTA.ArrElem[p.iPhi][p.iEta][p.iInvpt].nHitLayers;		
				H2nHitVSnMin.Fill(minLayers,nHitLayers);
				
				if(!par.gen_TrackFit)foundTrack = false;
				
				else {							
						double chi2, phi, eta, invPt, z0, t0, mass;
						
						if(verbose) cout << "fit candidate " << iCan << " [" << p.iPhi << "][" << p.iEta << "][" << p.iInvpt << "]" << endl;
						
						if(verbose) HTA.ArrElem[p.iPhi][p.iEta][p.iInvpt].printCandidate(ev, cout);		
				
						int nHits;				
						int retCodeFit = HTA.ArrElem[p.iPhi][p.iEta][p.iInvpt].
												fitCandidate(ev, chi2, phi, eta, invPt, z0, t0, mass, nLayers);
												
						if(verbose) cout << "nLayers: " << nLayers << endl;
												
						if(verbose) if(retCodeFit == -1) cout << " no fit" << endl;
						
						if(retCodeFit == 0) {
							HFitChi2.Fill(chi2);
							HFitLayers.Fill(nLayers);
							switch(nLayers) {
							  case 5 :
								HFit5Chi2.Fill(chi2);
								break;
							  case 6 :
								HFit6Chi2.Fill(chi2); 
								break;						  		
							  case 7 :
								HFit7Chi2.Fill(chi2); 
								break;
							  case 8 :
								HFit8Chi2.Fill(chi2);					 
						   }
						}
						
						if( (retCodeFit == 0) && (chi2 <= par.gen_chi2Cut) ) {//this is a good fit
																			
							foundTrack = true;	
							
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
				
							bool done =  HTA.ArrElem[i][j][k].thisCellDone;
							
							if(!done){ // this is executed only for one copy of duplicates
							
								HTA.ArrElem[i][j][k].thisCellDone = true;
							
								FitTrack ft;
									ft.chi2 = chi2;
									ft.phi = phi;
									ft.eta = eta;
									ft.invPt = invPt;
									ft.z0 = z0;
									ft.t0 = t0;
									ft.mass = mass;
									ft.nLayers = nLayers;
																						
								foundTracks.push_back(ft);
						
								if(verbose) cout << " good fit"	 << endl;					
								
							
								int DoF = 3*nLayers - 5;			
								if(verbose) cout << " DoF " << DoF << " chi2 " << chi2 << endl;
								if(verbose) cout << " result  pt:" << 1./invPt << " eta: " << eta << " phi: " << phi << " z0: " << z0 << " t0:" << t0 << endl;					
								if(ev.trackList.size()) 
									if(verbose) cout << "  track  pt:" << 1./ev.trackList[0].invPt << " eta: " << ev.trackList[0].eta 
									<< " phi: " << ev.trackList[0].phi << " z0: " << ev.trackList[0].z0 << " t0:" << ev.trackList[0].t0 << endl;
							}// end trick to eliminate duplicates
													
						}// end this is a good fit (retcode and chi2)
								
					}// end if TrackFit							
								
			} // end loop on candidates
			
			candidateRate += (double)nCan;			
						
			if(foundTrack) ++nEventsWithTracks;
			
			
	
		unsigned nFoundTracks = foundTracks.size();
		HnFoundTracks.Fill(nFoundTracks);
		
		if((nFoundTracks>0) && (ev.trackList.size()>0)){	
		
			// Fill efficiency histograms 				
		
			Track thisTrack = ev.trackList[0];
			HTrackEtaEff.Fill(thisTrack.eta);
			HTrackPhiEff.Fill(thisTrack.phi);
			HTrackInvPtEff.Fill(thisTrack.invPt);	
		}
		
		cout << nCan << " candidates and " << nFoundTracks << " tracks found in this event" << endl;
		if(ev.trackList.size()) 
				cout << "  track  pt:" << 1./ev.trackList[0].invPt << " eta: " << ev.trackList[0].eta 
					<< " phi: " << ev.trackList[0].phi << " z0: " << ev.trackList[0].z0 << " t0: " << ev.trackList[0].t0 << " mass: " << ev.trackList[0].mass << endl;
					
		for(unsigned iFT = 0; iFT != nFoundTracks; ++iFT)
				cout << " result  pt:" << 1./foundTracks[iFT].invPt << " eta: " << foundTracks[iFT].eta << " phi: " 
					<< foundTracks[iFT].phi << " z0: " << foundTracks[iFT].z0 << " t0:" << foundTracks[iFT].t0 
						<< " mass: " << foundTracks[iFT].mass  << " nLayers: " << foundTracks[iFT].nLayers << " chi2: " << foundTracks[iFT].chi2 << endl;			
								
		
		if(nFoundTracks == 0){
			if(ev.trackList.size()){
				unsigned i,j,k;
				ev.trackList[0].getIJK(i,j,k);
				cout << "cell [" << i << "][" << j << "][" << k << "] " ;
				cout << HTA.ArrElem[i][j][k].nHitLayers << " layers hit,";
				cout << " minLayers: " << HTA.ArrElem[i][j][k].minLayers << endl;
			}
		};
		
		cout << nEventsWithTracks << " event with tracks found out of " << iEv + 1 << " events" << endl;	
	
		///////////// END SUMMARY////////////////////////////////////////
			
		
		// Fill histogram of number of hits in each HTA cell 

		 for(int i = 0; i != HTA_NphiBins; ++i) 
				for(int j = 0; j != HTA_NetaBins; ++j)
					for(int k = 0; k != HTA_NinvptBins; ++k){										
						for(auto it = HTA.ArrElem[i][j][k].layerIndHitStat.begin(); it != HTA.ArrElem[i][j][k].layerIndHitStat.end();++it ){     						
								int n = (it->second).nHits;
								HCellStat.Fill(n);
								//cout << j << " " << k << " " << n << endl;					 						
						} 
					}

	
	/////////////////////////////////////////////////////////////////////////////////
	// Fill parameter resolution histograms
	
		if(ev.trackList.size()){
		
			double phi = ev.trackList[0].phi;
			double eta = ev.trackList[0].eta;
			double invPt = ev.trackList[0].invPt;
			double z0 = ev.trackList[0].z0;
			double t0 = ev.trackList[0].t0; 
		
			for( unsigned iT = 0; iT != nFoundTracks; ++iT){
		
				double DPhi = foundTracks[iT].phi - phi;
				double DEta = foundTracks[iT].eta - eta;
				double DPt  = foundTracks[iT].invPt - invPt;
				double DZ0  = foundTracks[iT].z0 - z0;
				double DT0  = foundTracks[iT].t0 - t0;
							
				double D2Phi = DPhi*DPhi;			
				double D2Eta = DEta*DEta;
				double D2Pt  = DPt*DPt;
											
				double D4Phi = D2Phi*D2Phi;
				double D4Eta = D2Eta*D2Eta;
				double D4Pt  = D2Pt*D2Pt;
				
				double P = fabs(sinh(eta)/invPt); // momentum
				if(P <= par.gen_massFitMaxP) 						
					HTrackMass.Fill(foundTracks[iT].mass);		
		
			   	HDeltaPhi.Fill(DPhi);
			   	HDeltaEta.Fill(DEta);
			   	HDeltaInvPt.Fill(DPt);
			   	HDeltaZ0.Fill(DZ0);
			   	HDeltaT0.Fill(DT0);	
			   	
			   	HPhiN.Fill(phi);
			   	HEtaN.Fill(eta);
			   	HInvPtN.Fill(invPt);
			   	
				HBiasPhiVsPhi.Fill(phi,DPhi);
				HBiasEtaVsPhi.Fill(phi,DEta);
				HBiasInvPtVsPhi.Fill(phi,DPt);

				HBiasPhiVsEta.Fill(eta,DPhi);
				HBiasEtaVsEta.Fill(eta,DEta);
				HBiasInvPtVsEta.Fill(eta,DPt);

				HBiasPhiVsInvPt.Fill(invPt,DPhi);
				HBiasEtaVsInvPt.Fill(invPt,DEta);
				HBiasInvPtVsInvPt.Fill(invPt,DPt);
				
				
				// HSigmaxxx histograms temporarily hold 
				// sum of delta squared of parameters
			
				HSigmaPhiVsPhi.Fill(phi,D2Phi);
				HSigmaEtaVsPhi.Fill(phi,D2Eta);
				HSigmaInvPtVsPhi.Fill(phi,D2Pt);

				HSigmaPhiVsEta.Fill(eta,D2Phi);
				HSigmaEtaVsEta.Fill(eta,D2Eta);
				HSigmaInvPtVsEta.Fill(eta,D2Pt);

				HSigmaPhiVsInvPt.Fill(invPt,D2Phi);
				HSigmaEtaVsInvPt.Fill(invPt,D2Eta);
				HSigmaInvPtVsInvPt.Fill(invPt,D2Pt);
				
				
				HD4PhiVsPhi.Fill(phi,D4Phi);
				HD4EtaVsPhi.Fill(phi,D4Eta);
				HD4InvPtVsPhi.Fill(phi,D4Pt);

				HD4PhiVsEta.Fill(eta,D4Phi);
				HD4EtaVsEta.Fill(eta,D4Eta);
				HD4InvPtVsEta.Fill(eta,D4Pt);

				HD4PhiVsInvPt.Fill(invPt,D4Phi);
				HD4EtaVsInvPt.Fill(invPt,D4Eta);
				HD4InvPtVsInvPt.Fill(invPt,D4Pt);
	
	
	
			}	
	
		}
	
	
	
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
	cout << "Number of candidates per Collision = " << candidateRate/nEvents*rateScale << endl;
	
	
	
	
	
	
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// Create resolution plots as a function of track parameters
{	
	
	TH1D * nevHist;
	TH1D * sumHist;
	TH1D * sum2Hist;
	TH1D * sum4Hist;
	
	for(int iCase = 0; iCase != 9; ++iCase){
			
		switch(iCase){
		
			case 0:
				nevHist  = &HPhiN;
				sumHist  = &HBiasPhiVsPhi;				
				sum2Hist = &HSigmaPhiVsPhi;	
				sum4Hist = &HD4PhiVsPhi;
				break;
				
			case 1:
				nevHist  = &HPhiN;
				sumHist  = &HBiasEtaVsPhi;
				sum2Hist = &HSigmaEtaVsPhi;
				sum4Hist = &HD4EtaVsPhi;
				break;
				
			case 2:
				nevHist  = &HPhiN;
				sumHist  = &HBiasInvPtVsPhi;
				sum2Hist = &HSigmaInvPtVsPhi;
				sum4Hist = &HD4InvPtVsPhi;
				break;
				
			case 3:
				nevHist  = &HEtaN;
				sumHist  = &HBiasPhiVsEta;
				sum2Hist = &HSigmaPhiVsEta;
				sum4Hist = &HD4PhiVsEta;
				break;
				
			case 4:
				nevHist  = &HEtaN;
				sumHist  = &HBiasEtaVsEta;
				sum2Hist = &HSigmaEtaVsEta;
				sum4Hist = &HD4EtaVsEta;
				break;
				
			case 5:
				nevHist  = &HEtaN;
				sumHist  = &HBiasInvPtVsEta;
				sum2Hist = &HSigmaInvPtVsEta;
				sum4Hist = &HD4InvPtVsEta;
				break;
				
			case 6:
				nevHist  = &HInvPtN;
				sumHist  = &HBiasPhiVsInvPt;
				sum2Hist = &HSigmaPhiVsInvPt;
				sum4Hist = &HD4PhiVsInvPt;
				break;
				
			case 7:
				nevHist  = &HInvPtN;
				sumHist  = &HBiasEtaVsInvPt;
				sum2Hist = &HSigmaEtaVsInvPt;
				sum4Hist = &HD4EtaVsInvPt;
				break;
						
			case 8:
				nevHist  = &HInvPtN;
				sumHist  = &HBiasInvPtVsInvPt;
				sum2Hist = &HSigmaInvPtVsInvPt;
				sum4Hist = &HD4InvPtVsInvPt;
				break;	
		
			}
			
			unsigned nBins0 = nevHist->GetNbinsX();	
			unsigned nBins1 = sumHist->GetNbinsX();	
			unsigned nBins2 = sum2Hist->GetNbinsX();	
			unsigned nBins4 = sum4Hist->GetNbinsX();
			
			if((nBins0 != nBins1)||(nBins0 != nBins2)||(nBins0 != nBins4))
				cout << "***** ERROR ***** nBins mismatch. Case " << iCase << " Efficiency plots skipped ";
			
			else for(int i = 1; i != nBins0+1; ++i){
			
				double sum = sumHist->GetBinContent(i);
				double nev =  nevHist->GetBinContent(i);
				double bias = sum/nev;
				
				double sumSquares = sum2Hist->GetBinContent(i);
				double sigma2 = sumSquares/nev - bias*bias;
				double sigma = sqrt(sigma2);
				
				double sumFourths = sum4Hist->GetBinContent(i);
				double sigma4 = sumFourths/nev;
				
				double sigmaErr = sqrt((sigma4 - sigma2*sigma2)/nev);
				 		sigmaErr =  sigmaErr/sigma/2.;
				/* 		
				 cout << sigma4 << "   ";		
				 cout << (sigma4 - sigma2*sigma2) << "   ";	
				  cout << (sigma4 - sigma2*sigma2)/sum << endl;
				 		
				cout << " iCase " << iCase << " " << i << " sigma " << sigma << " sigma2 " << sigma2 << " sigmaErr " << sigmaErr << endl;
				*/
				
				//double sigmaErr = sigma/sqrt(2*nev);			
				
				sumHist->SetBinContent (i, bias);	
				sumHist->SetBinError (i, sigma/sqrt(nev));
							
				sum2Hist->SetBinContent (i, sigma);	
				sum2Hist->SetBinError (i, sigmaErr);	
				
				sum4Hist->SetBinContent (i, sigma4);	
				sum4Hist->SetBinError (i, 0.);	
								
				 
			}
	
	}
}// end create resolution plots
	

/////////////////////////////////////////////////////////////////////////////////	
//
//  Create efficiency plots as a function of track parameters
//
{

	TGraph * g_eff1 = makeEffGraphFromHists(&HTrackPhiEff,&HTrackPhi);	
  	g_eff1->SetTitle("efficiency vs phi");
  	g_eff1->SetName("effVsPhi");  	
  	g_eff1->Write();
  	  
	TGraph * g_eff2 = makeEffGraphFromHists(&HTrackEtaEff,&HTrackEta);
  	g_eff2->SetTitle("efficiency vs eta");
  	g_eff2->SetName("effVsEta");
  	g_eff2->Write();

	TGraph * g_eff3 = makeEffGraphFromHists(&HTrackInvPtEff,&HTrackInvPt);
  	g_eff3->SetTitle("efficiency vs invPt");
  	g_eff3->SetName("effVsInvPt");
  	g_eff3->Write();
  	
} // end create efficiency plots 


/////////////////////////////////////////////////////////////////////////////////	
//
//  Create plots of the average number of hits per track as a function of eta and PT
//
{

// --------- Eta------------------------------------------------------------------- 
  {
	TH1D * nevHist;
	TH1D * sumHist;
	TH1D * sum2Hist;
	
	nevHist  = &HTrackEta;
	sumHist  = &HTrackNhitsVsEta;
	sum2Hist  = &HTrackNhits2VsEta;
		
	unsigned nBins0 = nevHist->GetNbinsX();	
	unsigned nBins1 = sumHist->GetNbinsX();		
	unsigned nBins2 = sum2Hist->GetNbinsX();						

	if((nBins0 != nBins1)||(nBins0 != nBins2))
		cout << "***** ERROR ***** nBins mismatch. Nhits vs. eta plot skipped ";
	
	else for(int i = 1; i != nBins0+1; ++i){		
		double nev =  nevHist->GetBinContent(i);
		if((unsigned)nev == 0) sumHist->SetBinContent(i,0.);
		else{
			double sum = sumHist->GetBinContent(i);
			double mean = sum/nev;
			sumHist->SetBinContent(i,mean);
			double sum2 = sum2Hist->GetBinContent(i);
			double err = sqrt(sum2/nev - mean*mean);
			sumHist->SetBinError(i,err);
		}
	}
  } // end number oh hits as a function of eta
  
// --------- PT------------------------------------------------------------------- 
	
  {
	TH1D * nevHist;
	TH1D * sumHist;
	TH1D * sum2Hist;
	
	nevHist  = &HTrackInvPt;
	sumHist  = &HTrackNhitsVsInvPt;
	sum2Hist  = &HTrackNhits2VsInvPt;
		
	unsigned nBins0 = nevHist->GetNbinsX();	
	unsigned nBins1 = sumHist->GetNbinsX();		
	unsigned nBins2 = sum2Hist->GetNbinsX();						

	if((nBins0 != nBins1)||(nBins0 != nBins2))
		cout << "***** ERROR ***** nBins mismatch. Nhits vs. PT plot skipped ";
	
	else for(int i = 1; i != nBins0+1; ++i){		
		double nev =  nevHist->GetBinContent(i);
		if((unsigned)nev == 0) sumHist->SetBinContent(i,0.);
		else{
			double sum = sumHist->GetBinContent(i);
			double mean = sum/nev;
			sumHist->SetBinContent(i,mean);
			double sum2 = sum2Hist->GetBinContent(i);
			double err = sqrt(sum2/nev - mean*mean);
			sumHist->SetBinError(i,err);
		}
	}
  } // end hits as a function of pt
  	

} // end plots of the average number of hits

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

		
	cout << "Writing histogram file..." << endl;
	histFile->Write(); // write histogram file
	
	
	
	return 0;		
	
} // end main



