//
//  Parameters.h
//  MuonColliderToy
//
//  Created by Luciano Ristori on 4/25/23
//


// dimensions of the HT array
// (global variables)

const static unsigned HTA_NphiBins = 15;
const static unsigned HTA_NetaBins = 360;
const static unsigned HTA_NinvptBins = 6;
		

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
	
	unsigned train_nTracksPerCell = 500; 
	
// Histogram files for three training phase

	string train_histFileName1 = "AAATrainingHists1.root";
	string train_histFileName2 = "AAATrainingHists2.root";
	string train_histFileName3 = "AAATrainingHists3.root";
	
// File containing all the data to init the Hough Transform array
// This is where all the data from training are written to
	
	string train_dataFileName = "HTAdata.txt";
	
	
// seed random generator

	long int train_randomSeed = 123456;
	long int train_randomSeed_trk = 0; // for track generation;

//////////////////////////////////////////////////////////////////////// 
// EVENT GENERATION SECTION  ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// 

// Draw plots on the screen at end of job

	bool gen_drawPlots = false;

// Number of events to be generated for simulation

	unsigned gen_nEvents = 100;

// number of tracks to be generated for each event

	unsigned gen_nTracks = 1; // Number of tracks per event
	
// fraction of BIB to be generated

	double gen_fracBib = -1.; // -1. == no BIB hit
	
// perform track fitting of candidates

	bool gen_TrackFit = true;
	unsigned gen_minLayersForFit = 5;
	double gen_chi2Cut = 45.;
	bool gen_massFit = false;
	double gen_massFitMaxP = 3.; // GeV/c
	
// limits in phi for BIB generation
	
	 //double gen_phia = -0.00216988; // for 500 mrad phi sector
	 //double gen_phib = 0.806467;    //
	 
	 //double gen_phia = -0.00312222; // for 100 mrad phi sector
 	 //double gen_phib = 0.406729;    //
 	  	 
 	//double gen_phia = -0.; // for 40 mrad phi sector
 	//double gen_phib = 0.346778;    // 	
 	
 	double gen_phia = -0.306759;// for 524 mrad phi sector
 	double gen_phib = 0.830405; // +- charges

 
 	 
// print candidates for every event

	bool gen_printCandidates = true;

// optimization mode for HTA fill = 0 (safe and slow) or 1 (faster) or 2 (fastest)
	
	int gen_fillMode = 2; 
	
// max error for eta projection in fillMode = 1 and 2

	int gen_errEta = 1;
	
// max error for phi projection in fillMode = 2

	int gen_errPhi = 1;
			
// create a file with data to 3D plot track candidates
// this file is the input for PlotTracks

	bool gen_PlotTracks = false;	
	string gen_plotDataFileName = "PlotData.txt";

// random generator seeds 

	long int gen_randomSeed = 121348;
	long int gen_randomSeed_trk = 220472; // for track generation
	
// data file from where HT array is initialized from

	string gen_dataFileName = "HTAdata_15_360_6_BigEtaSliceNew-Pt_3p0-LongBarrel.txt";
	//string gen_dataFileName = "HTAdata_15_360_6_BigEtaSliceNew.txt";
	//string gen_dataFileName = "HTAdata_15_360_6_BigEtaSliceNew-Pt_3p0-TimeRes_120ps_Vertex.txt";
	
	//string gen_dataFileName = "HTAdata.txt";
	
// data file with all Bib hits

	string gen_bibFileName = "BIBdata.txt";
	
// event generation histogram file

	string gen_histFileName = "AAAEventGeneration.root";
	//string gen_histFileName = "AAAEventGeneration_massFit_10ps.root";
	//string gen_histFileName = "AAA1nsecEventGeneration.root";
	//string gen_histFileName = "AAAAParticleIDTests.root";
	
	bool gen_fillHitHistograms = true; // selective histogram filling
	
// printing control

	bool gen_verbose = false;

//////////////////////////////////////////////////////////////////////// 
// MISCELLANEA /////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// 


// magnetic field in Tesla

	double magneticField = 4.0;

// coordinates of a single HTA cell and detector layer 
// whose hit coordinates we want to plot in 3-D

	int HTA_plotBinX = 5;
	int HTA_plotBinY = 5;
	int HTA_plotBinZ = 5;
	int HTA_plotLay = 5;


//////////////////////////////////////////////////////////////////////// 
// MEASUREMENT ERRORS //////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// 

	// detector measurement errors

	double geo_xphiSigmaB = 0.01;// Barrels: sigma x in mm
	double geo_xphiSigmaD = 0.01;// Discs: sigma x in mm
	double geo_zSigmaB = 0.1;// Barrels: sigma z in mm
	double geo_rSigmaD = 0.1;// Discs: sigma R in mm
	
	//double geo_tSigmaB = 300.;// Barrels: sigma t -  300 mm == 1000 ps
	//double geo_tSigmaD = 300;// Discs: sigma t - 300 mm == 1000 ps
		
	//double geo_tSigmaB = 36.;// Barrels: sigma t -  36 mm == 120 ps
	//double geo_tSigmaD = 36.;// Discs: sigma t - 36 mm == 120 ps

	double geo_tSigmaB = 18.;// Barrels: sigma t -  18 mm == 60 ps
	double geo_tSigmaD = 18.;// Discs: sigma t - 18 mm == 60 ps
	
	//double geo_tSigmaB = 3.;// Barrels: sigma t -  3 mm == 10 ps
	//double geo_tSigmaD = 3.;// Discs: sigma t - 3 mm == 10 ps
	
	//double geo_tSigmaB = 0.3;// Barrels: sigma t -  0.3 mm == 1 ps
	//double geo_tSigmaD = 0.3;// Discs: sigma t - 0.3 mm == 1 ps
	
	double geo_ineffB = 0.01;// Barrel inefficiency
	double geo_ineffD = 0.01;// Disc inefficiency


//////////////////////////////////////////////////////////////////////// 
// TRACK PARAMETERS ////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// 


	// track parameters - default values

	 double geo_def_t_phi = 0.262; // track phi mean
	 double geo_def_t_deltaPhi = 0.262; // track delta phi
	
	 double geo_def_t_eta = 0.; // track eta mean
	 double geo_def_t_deltaEta = 2.5; // track delta eta
	 double geo_def_t_invPt_max = 1./3.; // track invPt max Gev/c^(-1)  
	 double geo_def_t_invPt_min = -1./3; // track invPt min Gev/c^(-1)
	 double geo_def_t_x0 = 0.0; // track mean x0
	 double geo_def_t_y0 = 0.0; // track mean y0 
	 double geo_def_t_z0 = 0.0; // track mean z0
	 double geo_def_t_t0 = 0.0; // track mean t0
	 double geo_def_t_deltaX0 = 0.0; // track sigma x0
	 double geo_def_t_deltaY0 = 0.0; // track sigma y0
	 double geo_def_t_deltaZ0 = 1.5; // track sigma z0
	 double geo_def_t_deltaT0 = 1.5; // track sigma t0 mm
	
	
	// expansion factors applied to default parameter spaces
	// for EventGeneration(gen), HTA, and HTA training (train)
	
	 double gen_fiducial = 1.0;
	 double HTA_fiducial = 1.0;
	 double train_fiducial = 1.0;


	// track parameters for event generation

	 bool geo_gen_default = true; // apply default values?

	// if not use the following:
		// momentum range for mass plots
	 double geo_gen_t_phi = 0.262; // track phi mean
	 double geo_gen_t_deltaPhi = 0.262; // track delta phi
	 double geo_gen_t_eta = 0.; // track eta mean
	 double geo_gen_t_deltaEta = 2.0; // track delta eta
	 double geo_gen_t_invPt_max = 1./1.5; // track invPt max Gev/c^(-1)  
	 double geo_gen_t_invPt_min = 1./3.; // track invPt min Gev/c^(-1)
	 double geo_gen_t_x0 = 0.0; // track mean x0
	 double geo_gen_t_y0 = 0.0; // track mean y0 
	 double geo_gen_t_z0 = 0.0; // track mean z0
	 double geo_gen_t_t0 = 0.0; // track mean t0
	 double geo_gen_t_deltaX0 = 0.0; // track sigma x0
	 double geo_gen_t_deltaY0 = 0.0; // track sigma y0
	 double geo_gen_t_deltaZ0 = 1.5; // track sigma z0
	 double geo_gen_t_deltaT0 = 1.5; // track sigma t0 mm
	


	// track parameters for Hough Transform Array

	const bool HTA_default = true; // apply default values?

	// if not use the following:

	double HTA_t_phi = 0.005; // track phi mean
	double HTA_t_deltaPhi = 0.005; // track delta phi
	double HTA_t_eta = 1.0; // track eta mean
	double HTA_t_deltaEta = 1.0; // track delta eta
	double HTA_t_invPt_max = 1./3.; // track invPt max Gev/c^(-1)  
	double HTA_t_invPt_min = -1./3.; // track invPt min Gev/c^(-1)
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
	double geo_train_t_invPt_min = -1./3.; // track invPt min Gev/c^(-1)
	double geo_train_t_x0 = 0.0; // track mean x0
	double geo_train_t_y0 = 0.0; // track mean y0 
	double geo_train_t_z0 = 0.0; // track mean z0
	double geo_train_t_t0 = 0.0; // track mean t0
	double geo_train_t_deltaX0 = 0.0; // track sigma x0
	double geo_train_t_deltaY0 = 0.0; // track sigma y0
	double geo_train_t_deltaZ0 = 1.5; // track sigma z0
	double geo_train_t_deltaT0 = 1.5; // track sigma t0 mm
	
	


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
