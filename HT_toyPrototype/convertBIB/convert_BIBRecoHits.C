#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include <algorithm>
#include <random>
#include <fstream> 


#include "Geometry_Detector.cpp"


const float GeVtoMeV = 1000.;
const double c = 299.792458; // mm/ns - speed of light

TFile *input_file;
TTree *trkHits;
TFile* histFile = new TFile("AAAhists.root","RECREATE");  // histogram file
TH2I *H2HitRhoVsZa = new TH2I("HitRhoVsZa","HitRhoVsZa",5000,-2500.,+2500., 1800, 0., 1800.);


class Record {
	public:
	 	char hitType ; // B (barrel) or D (disc)
 	 	unsigned iLayer ; // indexed separately for barrels and discs
 	 	double x2 ; // secondary coordinate (mm)
 	 	double time; //  time (mm)
 	 
	void print(std::ostream &out){
		out << "type: " << hitType << " iLayer: " << iLayer << " x2: " << x2 << " time: " << time << endl;
	} 
	
	void write(std::ostream &out){
		out << hitType << " " << iLayer << " " << x2 << " " << time << endl;
	} 
};

vector<Record> buffer;
Record dummyRecord;


void convert_BIBRecoHits(const TString filename="allHits_ntuple_BIB.root"){

  Geometry g;
  g.print(cout);

  input_file = new TFile(filename, "read");
  
  trkHits  = (TTree*) input_file->Get("TrackerHitsTuple");

  
  // ================================================================================================
  // --- tracker hits

  int trk_n; 
  int *trk_id0 = new int[10000000]; 
  double *trk_x = new double[10000000]; 
  double *trk_y = new double[10000000];
  double *trk_z = new double[10000000];
  float *trk_t = new float[10000000];
  float *trk_e = new float[10000000];

  trkHits->SetBranchAddress("ntrh",&trk_n);// number of hits
  trkHits->SetBranchAddress("thci0",trk_id0);// hit code
  trkHits->SetBranchAddress("thpox",trk_x);// hit position x
  trkHits->SetBranchAddress("thpoy",trk_y);// hit position y
  trkHits->SetBranchAddress("thpoz",trk_z);// hit position z
  trkHits->SetBranchAddress("thtim",trk_t);// hit time
  trkHits->SetBranchAddress("thedp",trk_e);// hit deposited energyz	
  
  
  
  // hit parameters to be exported to toy MC
  char hitType = '0'; // B (barrel) or D (disc)
  unsigned iLayer = 0; // indexed separately for barrels and discs
  double x1 = 0., x2 = 0., time = 0.; // primary coordinate (phi), secondary (z or r), time (all in mm)
  
  
  unsigned maxTotHits = 10000000;
  //unsigned maxTotHits = 100;
  
  unsigned totHits = 0;
    
  for(int ientry=0; ientry<trkHits->GetEntries(); ++ientry){
	//for(int ientry=0; ientry<1; ++ientry){

	trkHits->GetEntry(ientry);

    for (int ihit=0; ihit<trk_n; ++ihit){
    	//cout << ihit << endl;
    
    	//if((ihit % 1000000) == 0) std::cout << ihit << " hits processed" << endl;

	  	// CellID encoding: "system:5,side:-2,layer:6,module:11,sensor:8"
		const unsigned int system = (unsigned) ( trk_id0[ihit] & 0x1f );
		const int side = (int) ( (trk_id0[ihit] >> 5) & 0x3 );
		const unsigned int layer = (unsigned) ( (trk_id0[ihit] >> 7) & 0x3f );
		const unsigned int module = (unsigned) ( (trk_id0[ihit] >> 13) & 0x7ff );
		const unsigned int sensor = (unsigned) ( (trk_id0[ihit] >> 24) & 0xff );

	  	float rho = sqrt(trk_x[ihit]*trk_x[ihit]+trk_y[ihit]*trk_y[ihit]); 
	  
	 	H2HitRhoVsZa->Fill(trk_z[ihit], rho);
	 	
	 	
	 	cout << endl << "s:"<< side << " z:" << trk_z[ihit] << " r:" << rho << " t:" << c*trk_t[ihit] << endl;
	 	
	 	if(side == 0){
	 		// Hit is in a barrel
	 		hitType = 'B';
	 		cout << "barrel " ;
	 		double minDelta = 1000.;
	 		int minB = -1;
	 		for(int iB = 0; iB != g.nBarrels; ++ iB){
	 			double delta = fabs(rho - g.B[iB].r);
	 			if(delta < minDelta) {
	 				minB = iB;
	 				minDelta = delta;
	 			}
	 		}
	 		iLayer = minB;
	 		cout << minB << " " << minDelta << endl;
	 		x2 = trk_z[ihit];
	 	
	 	}
	 	else if(side == 1){
	 		// Hit is in a right disc
	 		cout << "right disc ";
	 		hitType = 'D';
	 		double minDelta = 1000.;
	 		int nDhalf = g.nDiscs/2; // assume nDiscs is even
	 		int minD = -1;
	 		// right discs are from 0 to nDhalf - 1
	 		for(int iD = 0; iD != nDhalf; ++ iD){
	 			double delta = fabs(trk_z[ihit] - g.D[iD].z);
	 			if(delta < minDelta) {
	 				minD = iD;
	 				minDelta = delta;
	 			}
	 		}
	 		iLayer = minD;
	 		cout << minD << " " << minDelta << endl;
	 		x2 = rho;
	 	}
	 	else if(side == 3){
	 		// Hit is in a left disc
	 		cout << "left disc  " ;
	 		hitType = 'D';
	 		double minDelta = 1000.;
	 		int nDhalf = g.nDiscs/2; // assume nDiscs is even
	 		int minD = -1;
	 		// left discs are from nDhalf to nDiscs - 1
	 		for(int iD = nDhalf; iD != g.nDiscs; ++ iD){
	 			double delta = fabs(trk_z[ihit] - g.D[iD].z);
	 			if(delta < minDelta) {
	 				minD = iD;
	 				minDelta = delta;
	 			}
	 		}
	 		iLayer = minD;
	 		cout << minD << " " << minDelta << endl;
	 		x2 = rho;
	 	}
	 	
	 	x1 = trk_x[ihit];
	 	time = trk_t[ihit]*SpeedOfLight*1.e6;// conversion from ns to mm;
	 	
	 	cout << "tot: " << totHits << " type:"<< hitType << " layer:" << iLayer << " x1: " << x1 << " x2: " << x2 << " t: " << time << endl; 
	 	
	 	// check boundaries
	 	
	 	if(hitType == 'B'){   // Barrel: check z
	 		if(x2 < g.B[iLayer].zMin) {cout << "fail1" << endl; continue;}
	 		if(x2 > g.B[iLayer].zMax) {cout << "fail2" << endl; continue;}
	 	}
	 	else {              // Disc: check r
	 		if(x2 < g.D[iLayer].rMin) {cout << "fail3" << endl; continue;}
	 		if(x2 > g.D[iLayer].rMax) {cout << "fail4" << endl; continue;}	
	 	}
	 			
	 	dummyRecord.hitType = hitType;
	 	dummyRecord.iLayer = iLayer;
	 	dummyRecord.x2 = x2;
	 	dummyRecord.time = time;
	 	buffer.push_back(dummyRecord);
	 	
	 	
       	++totHits;
       	if(totHits >= maxTotHits) break;
     
    } // ihit loop
     
    if(totHits >= maxTotHits) break;
  
  } // trkHits ientry loop
  
  
   cout << endl;
   
   //for(unsigned i = 0; i != (unsigned)buffer.size(); ++i) buffer[i].print(cout);
   
   cout << "shuffling...." << endl;
   
   auto rng = std::default_random_engine {};
   std::shuffle(std::begin(buffer), std::end(buffer), rng);
   
  
  	cout << endl;
  	for(unsigned i = 0; i != (unsigned)buffer.size(); ++i) buffer[i].print(cout);
   
    // write data file
    	cout << endl;
		cout << "writing data file..." << endl;
		ofstream outfile;	
		outfile.open("BIBdata.txt");
		for(unsigned i = 0; i != (unsigned)buffer.size(); ++i) buffer[i].write(outfile);
		outfile.close();
   
	cout << endl;
	cout << buffer.size() << " records stored" << endl;
   

  delete [] trk_id0;
  delete [] trk_x;
  delete [] trk_y;
  delete [] trk_z;
  delete [] trk_t;
  delete [] trk_e;
    
    
  histFile->Write();
}
