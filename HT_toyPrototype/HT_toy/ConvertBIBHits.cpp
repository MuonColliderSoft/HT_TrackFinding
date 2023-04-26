#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include <algorithm>
#include <random>
#include <fstream> 

using namespace std;

#include "Geometry.cpp"


const float GeVtoMeV = 1000.;
const double c = 299.792458; // mm/ns - speed of light

const unsigned int verboseLevel = 0;

TFile *input_file;
TTree *trkHits;
TFile* histFile = new TFile("AAAConvertBIBHitsHists.root","RECREATE");  // histogram file
TH2I *H2HitRhoVsZa = new TH2I("HitRhoVsZa","HitRhoVsZa",5000,-2500.,+2500., 1800, 0., 1800.);


class Record {
public:
  char hitType ; // B (barrel) or D (disc)
  unsigned iLayer ; // indexed separately for barrels and discs
  double x2 ; // secondary coordinate (mm)
  double time; //  time (mm)
 	 
  void print(std::ostream &out){
    out << "type: " << hitType << " iLayer: " << iLayer << " x2: " << x2 << " time: " << time << std::endl;
  } 
	
  void write(std::ostream &out){
    out << hitType << " " << iLayer << " " << x2 << " " << time << std::endl;
  } 
};

std::vector<Record> buffer;
Record dummyRecord;


int main(int argc, char **argv){

  const char* rootfilename = ( argc > 1 ? argv[1] : "allHits_ntuple_BIB.root" );
  const char* outfilename = ( argc > 2 ? argv[2] : "BIBdata.txt");
  
  Geometry g;
  g.print(std::cout);

  input_file = new TFile(rootfilename, "read");
  
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
  
  unsigned totHits = 0;
    
  for(int ientry=0; ientry<trkHits->GetEntries(); ++ientry){

    trkHits->GetEntry(ientry);

    for (int ihit=0; ihit<trk_n; ++ihit){

      // CellID encoding: "system:5,side:-2,layer:6,module:11,sensor:8"
      const unsigned int system = (unsigned) ( trk_id0[ihit] & 0x1f );
      const int side = (int) ( (trk_id0[ihit] >> 5) & 0x3 );
      const unsigned int layer = (unsigned) ( (trk_id0[ihit] >> 7) & 0x3f );
      const unsigned int module = (unsigned) ( (trk_id0[ihit] >> 13) & 0x7ff );
      const unsigned int sensor = (unsigned) ( (trk_id0[ihit] >> 24) & 0xff );

      float rho = sqrt(trk_x[ihit]*trk_x[ihit]+trk_y[ihit]*trk_y[ihit]); 


      H2HitRhoVsZa->Fill(trk_z[ihit], rho);
      
	 	
      std::cout << std::endl << "s:"<< side << " z:" << trk_z[ihit] << " r:" << rho << " t:" << c*trk_t[ihit] << std::endl;
	 	
      if(side == 0){
	// Hit is in a barrel
	hitType = 'B';
	std::cout << "barrel " ;
	double minDelta = 1000.;
	int minB = -1;
	for(int iB = 0; iB != g.nBarrels; ++ iB){
	  if ( trk_z[ihit] < g.B[iB].zMin || trk_z[ihit] > g.B[iB].zMax ) continue;  
	  
	  double delta = fabs(rho - g.B[iB].r);
	  if(delta < minDelta) {
	    minB = iB;
	    minDelta = delta;
	  }
	}
	if ( minB == -1 ) continue;
	iLayer = minB;
	std::cout << minB << " " << minDelta << std::endl;
	x2 = trk_z[ihit];
	 	
      }
      else if(side == 1){
	// Hit is in a right disc
	std::cout << "right disc ";
	hitType = 'D';
	double minDelta = 1000.;
	int nDhalf = g.nDiscs/2; // assume nDiscs is even
	int minD = -1;
	// right discs are from 0 to nDhalf - 1
	for(int iD = 0; iD != nDhalf; ++ iD){
	  if ( rho < g.D[iD].rMin || rho > g.D[iD].rMax ) continue; 
	  double delta = fabs(trk_z[ihit] - g.D[iD].z);
	  if(delta < minDelta) {
	    minD = iD;
	    minDelta = delta;
	  }
	}
	if (  minD == -1 ) continue;
	iLayer = minD;
	std::cout << minD << " " << minDelta << std::endl;
	x2 = rho;
      }
      else if(side == 3){
	// Hit is in a left disc
	std::cout << "left disc  " ;
	hitType = 'D';
	double minDelta = 1000.;
	int nDhalf = g.nDiscs/2; // assume nDiscs is even
	int minD = -1;
	// left discs are from nDhalf to nDiscs - 1
	for(int iD = nDhalf; iD != g.nDiscs; ++ iD){
	  if ( rho < g.D[iD].rMin || rho > g.D[iD].rMax ) continue; 
	  double delta = fabs(trk_z[ihit] - g.D[iD].z);
	  if(delta < minDelta) {
	    minD = iD;
	    minDelta = delta;
	  }
	}
	if (  minD == -1 ) continue;
	iLayer = minD;
	std::cout << minD << " " << minDelta << std::endl;
	x2 = rho;
      }
	 	
      x1 = trk_x[ihit];
      time = trk_t[ihit]*SpeedOfLight*1.e6;// conversion from ns to mm;
	 	
      std::cout << "tot: " << totHits << " type:"<< hitType << " layer:" << iLayer << " x1: " << x1 << " x2: " << x2 << " t: " << time << std::endl; 
	 	
      // check boundaries
	 	
      if(hitType == 'B'){   // Barrel: check z
	if(x2 < g.B[iLayer].zMin) {std::cout << "fail1" << std::endl; continue;}
	if(x2 > g.B[iLayer].zMax) {std::cout << "fail2" << std::endl; continue;}
      }
      else {              // Disc: check r
	if(x2 < g.D[iLayer].rMin) {std::cout << "fail3" << std::endl; continue;}
	if(x2 > g.D[iLayer].rMax) {std::cout << "fail4" << std::endl; continue;}	
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
  
  
  std::cout << std::endl;
   
  std::cout << "shuffling...." << std::endl;
   
  auto rng = std::default_random_engine {};
  std::shuffle(std::begin(buffer), std::end(buffer), rng);
   
  
  std::cout << std::endl;
  if ( verboseLevel > 1 )
    for(unsigned i = 0; i != (unsigned)buffer.size(); ++i) buffer[i].print(std::cout);
   
  // write data file
  std::cout << std::endl;
  std::cout << "writing data file..." << std::endl;
  std::ofstream outfile;	
  outfile.open(outfilename);
  for(unsigned i = 0; i != (unsigned)buffer.size(); ++i) buffer[i].write(outfile);
  outfile.close();
   
  std::cout << std::endl;
  std::cout << buffer.size() << " records stored" << std::endl;
   

  delete [] trk_id0;
  delete [] trk_x;
  delete [] trk_y;
  delete [] trk_z;
  delete [] trk_t;
  delete [] trk_e;
    
    
  histFile->Write();

  return 0;

}
