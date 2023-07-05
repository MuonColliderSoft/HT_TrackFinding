//
//  Parameters.h
//  MuonColliderToy
//
//  Created by Luciano Ristori on 4/25/23
//


// dimensions of the HT array
// (global variables)

const static unsigned HTA_NphiBins = 1;
const static unsigned HTA_NetaBins = 200;
const static unsigned HTA_NinvptBins = 10;
		

class Parameters {

////////////////////////////////////////////////////////////////////////       
// TRAINING SECTION ////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// 
 
public:

// Hit coordinates distributions are diagonalized to find principal components
// before finding minimum and maximum  

	bool train_Diagonalize = false; 
	
// At the end of training, a summary of results is printed

	bool train_Summary = false;	
	
// number of events to be generated for training
	
	unsigned train_nEvents = 1e6; 
	
// Histogram files for three training phases

	string train_histFileName1 = "AAATrainingHists1.root";
	string train_histFileName2 = "AAATrainingHists2.root";
	string train_histFileName3 = "AAATrainingHists3.root";
	
// File containing all the data to init the Hough Transform array
// This is where all the data from training are written to
	
	string train_dataFileName = "HTAdata.txt";
	
// number of Poisson sigmas to accept as a down fluctuation of the number of hits
// in each detector layer to be accepted for a HTA cell
	
	double train_nSigmaCell = 3.;
	
// seed random generator

	long int train_randomSeed = 0;

//////////////////////////////////////////////////////////////////////// 
// EVENT GENERATION SECTION  ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// 

// Number of events to be generated for simulation

	unsigned gen_nEvents = 100;

// number of tracks to be generated for each event

	unsigned gen_nTracks = 0.; // Number of tracks per event
	
// fraction of BIB to be generated

	double gen_fracBib = 1.; // -1. == no BIB hit
	
// limits in phi for BIB generation

	double gen_phia = 0.0;
	double gen_phib = 0.344;
	
// optimization mode for HTA fill = 0 (safe and slow) or 1 (faster)	
	
	int gen_fillMode = 1; 
	
// max error for eta projection in fillMode = 1

	int gen_errEta = 2;
	
// create a file with data to 3D plot track candidates
// this file is the input for PlotTracks

	bool gen_PlotTracks = true;	
	string gen_plotDataFileName = "PlotData.txt";

// seed random generator

	long int gen_randomSeed = 121348;
	
// data file from where HT array is initialized from

	string gen_dataFileName = "HTAdata.txt";
	
// data file with all Bib hits

	string gen_bibFileName = "BIBdata.txt";
	
// event generation histogram file

	string gen_histFileName = "AAAEventGeneration.root";

//////////////////////////////////////////////////////////////////////// 
// MISCELLANEA /////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// 


// magnetic field in Tesla

	double magneticField = 4.0;

// coordinates of a single HTA cell and detector layer 
// whose hit coordinates we want to plot in 3-D

	int HTA_plotBinX = 0;
	int HTA_plotBinY = 100;
	int HTA_plotBinZ = 5;
	int HTA_plotLay = 7;


//////////////////////////////////////////////////////////////////////// 
// MEASUREMENT ERRORS //////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// 

	// detector measurement errors

	double geo_xphiSigmaB = 0.1;// Barrels: sigma x in mm
	double geo_xphiSigmaD = 0.1;// Discs: sigma x in mm
	double geo_zSigmaB = 0.1;// Barrels: sigma z in mm
	double geo_rSigmaD = 0.1;// Discs: sigma R in mm
	double geo_tSigmaB = 15.0;// Barrels: sigma t -  15 mm == 50 ps
	double geo_tSigmaD = 15.0;// Discs: sigma t - 15 mm == 50 ps


//////////////////////////////////////////////////////////////////////// 
// TRACK PARAMETERS ////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// 


	// track parameters - default values

	const double geo_def_t_phi = 0.005; // track phi mean
	const double geo_def_t_deltaPhi = 0.005; // track delta phi
	const double geo_def_t_eta = 1.0; // track eta mean
	const double geo_def_t_deltaEta = 1.0; // track delta eta
	const double geo_def_t_invPt_max = 1./3.; // track invPt max Gev/c^(-1)  
	const double geo_def_t_invPt_min = 0.; // track invPt min Gev/c^(-1)
	const double geo_def_t_x0 = 0.0; // track mean x0
	const double geo_def_t_y0 = 0.0; // track mean y0 
	const double geo_def_t_z0 = 0.0; // track mean z0
	const double geo_def_t_t0 = 0.0; // track mean t0
	const double geo_def_t_deltaX0 = 0.0; // track sigma x0
	const double geo_def_t_deltaY0 = 0.0; // track sigma y0
	const double geo_def_t_deltaZ0 = 1.5; // track sigma z0
	const double geo_def_t_deltaT0 = 1.5; // track sigma t0 mm
	
	
	// expansion factors applied to default parameter spaces
	// for EventGeneration(gen), HTA, and HTAtraining (train)
	
	const double gen_fiducial = 1.0;
	const double HTA_fiducial = 1.0;
	const double train_fiducial = 1.0;


	// track parameters for event generation

	const bool geo_gen_default = true; // apply default values?

	// if not use the following:

	double geo_gen_t_phi = 0.005; // track phi mean
	double geo_gen_t_deltaPhi = 0.005; // track delta phi
	double geo_gen_t_eta = 1.0; // track eta mean
	double geo_gen_t_deltaEta = 1.0; // track delta eta
	double geo_gen_t_invPt_max = 1./3.; // track invPt max Gev/c^(-1)  
	double geo_gen_t_invPt_min = 0.; // track invPt min Gev/c^(-1)
	double geo_gen_t_x0 = 0.0; // track mean x0
	double geo_gen_t_y0 = 0.0; // track mean y0 
	double geo_gen_t_z0 = 0.0; // track mean z0
	double geo_gen_t_t0 = 0.0; // track mean t0
	double geo_gen_t_deltaX0 = 0.0; // track sigma x0
	double geo_gen_t_deltaY0 = 0.0; // track sigma y0
	double geo_gen_t_deltaZ0 = 0.0; // track sigma z0
	double geo_gen_t_deltaT0 = 0.0; // track sigma t0 mm


	// track parameters for Hough Transform Array

	const bool HTA_default = true; // apply default values?

	// if not use the following:

	double HTA_t_phi = 0.005; // track phi mean
	double HTA_t_deltaPhi = 0.005; // track delta phi
	double HTA_t_eta = 1.0; // track eta mean
	double HTA_t_deltaEta = 1.0; // track delta eta
	double HTA_t_invPt_max = 1./3.; // track invPt max Gev/c^(-1)  
	double HTA_t_invPt_min = 0.; // track invPt min Gev/c^(-1)
	double HTA_t_x0 = 0.0; // track mean x0
	double HTA_t_y0 = 0.0; // track mean y0 
	double HTA_t_z0 = 0.0; // track mean z0
	double HTA_t_t0 = 0.0; // track mean t0
	double HTA_t_deltaX0 = 0.0; // track sigma x0
	double HTA_t_deltaY0 = 0.0; // track sigma y0
	double HTA_t_deltaZ0 = 0.0; // track sigma z0
	double HTA_t_deltaT0 = 0.0; // track sigma t0 mm


	// track parameters for training

	const bool geo_train_default = true; // apply default values?

	// if not use the following:

	double geo_train_t_phi = 0.005; // track phi mean
	double geo_train_t_deltaPhi = 0.005; // track delta phi
	double geo_train_t_eta = 1.0; // track eta mean
	double geo_train_t_deltaEta = 1.0; // track delta eta
	double geo_train_t_invPt_max = 1./3.; // track invPt max Gev/c^(-1)  
	double geo_train_t_invPt_min = 0.; // track invPt min Gev/c^(-1)
	double geo_train_t_x0 = 0.0; // track mean x0
	double geo_train_t_y0 = 0.0; // track mean y0 
	double geo_train_t_z0 = 0.0; // track mean z0
	double geo_train_t_t0 = 0.0; // track mean t0
	double geo_train_t_deltaX0 = 0.0; // track sigma x0
	double geo_train_t_deltaY0 = 0.0; // track sigma y0
	double geo_train_t_deltaZ0 = 0.0; // track sigma z0
	double geo_train_t_deltaT0 = 0.0; // track sigma t0 mm
	
	


//////////////////////////////////////////////////////////////////////// 
// CONSTRUCTOR//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// 

	Parameters(){
	
	// copy default values where appropriate
	
		if(geo_gen_default){
		
			geo_gen_t_phi = geo_def_t_phi; // track phi mean
			geo_gen_t_deltaPhi = geo_def_t_deltaPhi * gen_fiducial; // track delta phi
			geo_gen_t_eta = geo_def_t_eta; // track eta mean
			geo_gen_t_deltaEta = geo_def_t_deltaEta * gen_fiducial; // track delta eta
			geo_gen_t_invPt_max = geo_def_t_invPt_max * gen_fiducial; // track invPt max Gev/c^(-1)  
			geo_gen_t_invPt_min = geo_def_t_invPt_min; // track invPt min Gev/c^(-1)
			geo_gen_t_x0 = geo_def_t_x0; // track mean x0
			geo_gen_t_y0 = geo_def_t_y0; // track mean y0 
			geo_gen_t_z0 = geo_def_t_z0; // track mean z0
			geo_gen_t_t0 = geo_def_t_t0; // track mean t0
			geo_gen_t_deltaX0 = geo_def_t_deltaX0 * gen_fiducial; // track sigma x0
			geo_gen_t_deltaY0 = geo_def_t_deltaY0 * gen_fiducial; // track sigma y0
			geo_gen_t_deltaZ0 = geo_def_t_deltaZ0 * gen_fiducial; // track sigma z0
			geo_gen_t_deltaT0 = geo_def_t_deltaT0 * gen_fiducial; // track sigma t0 mm			
		
		}
	
		if(HTA_default){
		
			HTA_t_phi = geo_def_t_phi; // track phi mean
			HTA_t_deltaPhi = geo_def_t_deltaPhi * HTA_fiducial; // track delta phi
			HTA_t_eta = geo_def_t_eta; // track eta mean
			HTA_t_deltaEta = geo_def_t_deltaEta * HTA_fiducial; // track delta eta
			HTA_t_invPt_max = geo_def_t_invPt_max * HTA_fiducial; // track invPt max Gev/c^(-1)  
			HTA_t_invPt_min = geo_def_t_invPt_min; // track invPt min Gev/c^(-1)
			HTA_t_x0 = geo_def_t_x0; // track mean x0
			HTA_t_y0 = geo_def_t_y0; // track mean y0 
			HTA_t_z0 = geo_def_t_z0; // track mean z0
			HTA_t_t0 = geo_def_t_t0; // track mean t0
			HTA_t_deltaX0 = geo_def_t_deltaX0 * HTA_fiducial; // track sigma x0
			HTA_t_deltaY0 = geo_def_t_deltaY0 * HTA_fiducial; // track sigma y0
			HTA_t_deltaZ0 = geo_def_t_deltaZ0 * HTA_fiducial; // track sigma z0
			HTA_t_deltaT0 = geo_def_t_deltaT0 * HTA_fiducial; // track sigma t0 mm			
		
		}
	
	
		if(geo_train_default){
		
			geo_train_t_phi = geo_def_t_phi; // track phi mean
			geo_train_t_deltaPhi = geo_def_t_deltaPhi * train_fiducial; // track delta phi
			geo_train_t_eta = geo_def_t_eta; // track eta mean
			geo_train_t_deltaEta = geo_def_t_deltaEta * train_fiducial; // track delta eta
			geo_train_t_invPt_max = geo_def_t_invPt_max * train_fiducial; // track invPt max Gev/c^(-1)  
			geo_train_t_invPt_min = geo_def_t_invPt_min; // track invPt min Gev/c^(-1)
			geo_train_t_x0 = geo_def_t_x0; // track mean x0
			geo_train_t_y0 = geo_def_t_y0; // track mean y0 
			geo_train_t_z0 = geo_def_t_z0; // track mean z0
			geo_train_t_t0 = geo_def_t_t0; // track mean t0
			geo_train_t_deltaX0 = geo_def_t_deltaX0 * train_fiducial; // track sigma x0
			geo_train_t_deltaY0 = geo_def_t_deltaY0 * train_fiducial; // track sigma y0
			geo_train_t_deltaZ0 = geo_def_t_deltaZ0 * train_fiducial; // track sigma z0
			geo_train_t_deltaT0 = geo_def_t_deltaT0 * train_fiducial; // track sigma t0 mm			
		
		}
	}



};
