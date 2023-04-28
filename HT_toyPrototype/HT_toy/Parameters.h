//
//  Parameters.h
//  MuonColliderToy
//
//  Created by Luciano Ristori on 4/25/23
//


// dimensions of the HT array
//(global)

const static unsigned HTA_NphiBins = 1;
const static unsigned HTA_NetaBins = 200;
const static unsigned HTA_NinvptBins = 10;
		

class Parameters {

////////////////////////////////////////////////////////////////////////       
// HTArrayTraining.cpp /////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// 
 
public:

// Hit coordinates distributions are diagonalized to find principal components
// before finding minimum and maximum  

	bool train_Diagonalize = false; 
	
// At the end of training, a summary of results is printed

	bool train_Summary = false;	
	
// number of events to be generated	for training
	
	unsigned train_nEvents = 1e6; 
	
// Histogram files for three training phases

	string train_histFileName1 = "AAATrainingHists1.root";
	string train_histFileName2 = "AAATrainingHists2.root";
	string train_histFileName3 = "AAATrainingHists3.root";
	
// File containing all the data to init the Hough Transform array
// This is where all the data from training are written to
	
	string train_dataFileName = "HTAdata.txt";

//////////////////////////////////////////////////////////////////////// 
// EventGeneration.cpp  ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// 

// Number of events to be generated for simulation

	unsigned gen_nEvents = 100; 

// number of tracks to be generated for each event

	unsigned gen_nTracks = 1; // Number of tracks per event
	
// fraction of BIB to be generated

	double gen_fracBib = 0.;
	
// optimization mode for HTA fill = 0 (safe and slow) or 1 (fast)	
	
	int gen_fillMode = 1; 
	
// these should not be used by event generation	
	//bool Diagonalize = false;
	//bool Summary = false;
	
// create a file with data to 3D plot track candidates

	bool gen_PlotTracks = true;	
	string gen_plotDataFileName = "PlotData.txt";

// seed random generator

	long int gen_randomSeed = 121348;
	//generator.seed(121348); // seeding the random generator to be different from training
	
// data file from where HT array is initialized from

	string gen_dataFileName = "HTAdata.txt";
	
// data file with all Bib events

	string gen_bibFileName = "BIBdata.txt";	
	//int code = bibRead.readFile("BIBdata.txt");

// event generation histogram file

	string gen_histFileName = "AAAEventGeneration.root";
	//TFile* histFile = new TFile("AAAEventGeneration.root","RECREATE");  // histogram file


//////////////////////////////////////////////////////////////////////// 
// HTArray.cpp  ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// 


// coordinates of a single HTA cell and detector layer 
// whose hit coordinates we want to plot in 3-D

	int HTA_plotBinX = 0;
	int HTA_plotBinY = 100;
	int HTA_plotBinZ = 5;
	int HTA_plotLay = 7;

// this file holds information about how good is the approximation that
// we use to predict which cells can be actually relevant for a given hit
// and which ones we can safely skip. This file can be read by PlotMaps.cpp
// to create two meaningful histograms

	string HTA_mapDataFileName = "mapData.txt";
	//outfile.open("mapData.txt");


//////////////////////////////////////////////////////////////////////// 
// Geometry.cpp ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// 


// detector measurement errors

double geo_xphiSigmaB = 0.1;// Barrels: sigma x in mm
double geo_xphiSigmaD = 0.1;// Discs: sigma x in mm
double geo_zSigmaB = 0.1;// Barrels: sigma z in mm
double geo_rSigmaD = 0.1;// Discs: sigma R in mm
double geo_tSigmaB = 15.0;// Barrels: sigma t -  15 mm == 50 ps
double geo_tSigmaD = 15.0;// Discs: sigma t - 15 mm == 50 ps


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
const double geo_def_t_deltaZ0 = 0.0; // track sigma z0
const double geo_def_t_deltaT0 = 0.0; // track sigma t0 mm


// track parameters - training

const bool geo_train_default = true; // apply default values?

// if not use the following:

const double geo_train_t_phi = 0.005; // track phi mean
const double geo_train_t_deltaPhi = 0.005; // track delta phi
const double geo_train_t_eta = 1.0; // track eta mean
const double geo_train_t_deltaEta = 1.0; // track delta eta
const double geo_train_t_invPt_max = 1./3.; // track invPt max Gev/c^(-1)  
const double geo_train_t_invPt_min = 0.; // track invPt min Gev/c^(-1)
const double geo_train_t_x0 = 0.0; // track mean x0
const double geo_train_t_y0 = 0.0; // track mean y0 
const double geo_train_t_z0 = 0.0; // track mean z0
const double geo_train_t_t0 = 0.0; // track mean t0
const double geo_train_t_deltaX0 = 0.0; // track sigma x0
const double geo_train_t_deltaY0 = 0.0; // track sigma y0
const double geo_train_t_deltaZ0 = 0.0; // track sigma z0
const double geo_train_t_deltaT0 = 0.0; // track sigma t0 mm


// track parameters - event generation

const bool geo_gen_default = true; // apply default values?

// if not use the following:

const double geo_gen_t_phi = 0.005; // track phi mean
const double geo_gen_t_deltaPhi = 0.005; // track delta phi
const double geo_gen_t_eta = 1.0; // track eta mean
const double geo_gen_t_deltaEta = 1.0; // track delta eta
const double geo_gen_t_invPt_max = 1./3.; // track invPt max Gev/c^(-1)  
const double geo_gen_t_invPt_min = 0.; // track invPt min Gev/c^(-1)
const double geo_gen_t_x0 = 0.0; // track mean x0
const double geo_gen_t_y0 = 0.0; // track mean y0 
const double geo_gen_t_z0 = 0.0; // track mean z0
const double geo_gen_t_t0 = 0.0; // track mean t0
const double geo_gen_t_deltaX0 = 0.0; // track sigma x0
const double geo_gen_t_deltaY0 = 0.0; // track sigma y0
const double geo_gen_t_deltaZ0 = 0.0; // track sigma z0
const double geo_gen_t_deltaT0 = 0.0; // track sigma t0 mm


// track parameters - Hough Transform Array

const bool geo_HTA_default = true; // apply default values?

// if not use the following:

const double geo_HTA_t_phi = 0.005; // track phi mean
const double geo_HTA_t_deltaPhi = 0.005; // track delta phi
const double geo_HTA_t_eta = 1.0; // track eta mean
const double geo_HTA_t_deltaEta = 1.0; // track delta eta
const double geo_HTA_t_invPt_max = 1./3.; // track invPt max Gev/c^(-1)  
const double geo_HTA_t_invPt_min = 0.; // track invPt min Gev/c^(-1)
const double geo_HTA_t_x0 = 0.0; // track mean x0
const double geo_HTA_t_y0 = 0.0; // track mean y0 
const double geo_HTA_t_z0 = 0.0; // track mean z0
const double geo_HTA_t_t0 = 0.0; // track mean t0
const double geo_HTA_t_deltaX0 = 0.0; // track sigma x0
const double geo_HTA_t_deltaY0 = 0.0; // track sigma y0
const double geo_HTA_t_deltaZ0 = 0.0; // track sigma z0
const double geo_HTA_t_deltaT0 = 0.0; // track sigma t0 mm



};
