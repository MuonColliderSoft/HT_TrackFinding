#include <iostream> 
#include <fstream> 
#include <iomanip> 
#include <string>
#include <sstream> 
#include <algorithm> 
#include <random>
#include <list>

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"



#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"

#include "Statistics.cpp"



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

// random generator

//std::default_random_engine generator;
std::mt19937 generator; // Mersenne Twister 
std::mt19937 generator_trk; // Mersenne Twister for track generation
std::uniform_real_distribution<double> distribution(0.,1.);
std::normal_distribution<double> gauss(0.0,1.0);

//#include "Gauss.h"

#include "GlobalConstants.cpp"

#include "Parameters.h"
Parameters par; // class that holds all configurable parameters					   

#include "DetectorGeometry.cpp"
DetectorGeometry dg;

#include "TrackGeometry.cpp"
TrackGeometry tg("train"); // track parameters for HTA training

#include "Hit.h"
#include "Track.cpp"
#include "Hit.cpp"                                        
#include "BibFileReader.cpp"
#include "Event.cpp"
#include "HTArray.cpp"


bool Debug = false;
bool verbose = false;

int TrainingPhase = 0; // 0 = invalid 
		       // 1 = from scratch 
		       // 2 = second step after diagonalization
		       // 3 = verification and statistics
					   
					   
// The following global parameters are initialized in the main function
					   
unsigned nTracks; // number of tracks to be generated per cell for training		   
bool Diagonalize;
bool Summary;



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



int main(){
  

	// Parameter initialization:

	nTracks = par.train_nTracksPerCell;
	Diagonalize = par.train_Diagonalize;
	Summary = par.train_Summary;
	
	// seeding the random generators
  	long int randomSeed = par.train_randomSeed;
  	if(randomSeed)generator.seed(randomSeed); 
  	long int randomSeed_trk = par.train_randomSeed_trk;
  	if(randomSeed_trk)generator_trk.seed(randomSeed_trk); 
	
	
	// Instantiate Hough Transform Array
	
	cout << "Instantiating Hough Transform Array..." << endl;
	
	HTArray HTA;
	
	//cout << "Done" << endl;
	
	HTA.print(cout);
			
	// Open data file
	
	string dataFileName = par.train_dataFileName;
	cout << "Opening Data File: " << dataFileName << " ..." << endl;
	ifstream infile;
	infile.open(dataFileName);
	if(infile){ // data file found
		int retcode = HTA.read(infile);
		cout << "HTA.read retcode: " << retcode << endl;
		TrainingPhase = retcode;		
	} else {
		cout << "Data File not found, starting training from scratch" << endl;
		TrainingPhase = 1;
	}
	
	
	if(TrainingPhase < 0) return TrainingPhase; // datafile does not match current configuration

	string train_histFileName;
	switch (TrainingPhase) {
	case 1:
	  train_histFileName = par.train_histFileName1;
	  break;
	case 2:
	  train_histFileName = par.train_histFileName2;
	  break;
	case 3:
	  train_histFileName = par.train_histFileName3;
	  break;
	}

	TFile* histFile = new TFile(train_histFileName.c_str(),"RECREATE");  // histogram file
	
	HTA.initHists();
  
	dg.print(cout);
	
	
	TH1D HBarrelFraction("HBarrelFraction","HBarrelFraction", 11, -0.05, 1.05);HBarrelFraction.SetStats(false);
	TH2D HNDvsNB("NDvsNB","NDvsNB",  13, -0.5,+12.5, 13, -0.5, +12.5); HNDvsNB.SetStats(false);
	
	TH1D HCellStat("HCellStat","HCellStat",1000, 0., 5000.);
	TH1D Hd0Width("Hd0Width","Hd0Width",100,0.,100.);
	TH1D Hd1Width("Hd1Width","Hd1Width",100,0.,200.);
	TH1D Hd2Width("Hd2Width","Hd2Width",100,0.,200.);
	
	
	TH1D HitX("HitX","HitX", 4000, -2000.,+2000.); HitX.SetStats(false);
	TH1D HitY("HitY","HitY", 4000, -2000.,+2000.); HitY.SetStats(false);
	TH1D HitZ("HitZ","HitZ", 5000, -2500.,+2500.); HitZ.SetStats(false);
	TH2D HitXYbarrel("HitXYbarrel","HitXYbarrel", 4000, -2000.,+2000., 4000, -2000, +2000); HitXYbarrel.SetStats(false);
	TH2D HitXYdisc("HitXYdisc","HitXYdisc", 4000, -2000.,+2000., 4000, -2000, +2000); HitXYdisc.SetStats(false);
	TH2D HitRZ("HitRZ","HitRZ",  5000, -2500.,+2500., 1800, 0, 1800); HitRZ.SetStats(false);
	TH2I HitZPhi("HitZPhi","HitZPhi",  10000, -2500.,+2500., 10000, -1., +1.); HitZPhi.SetStats(false);
	TH2I HitRPhi("HitRPhi","HitRPhi",  10000, 0.,+2500., 10000, -1., +1.); HitRPhi.SetStats(false);
	
	
	TH1D HTrackZ0("HTrackZ0","HTrackZ0", 600, -300,+300); HTrackZ0.SetStats(true);
	TH1D HTrackT0("HTrackT0","HTrackT0", 600, -300,+300); HTrackT0.SetStats(true);
	TH1D HTrackEta("HTrackEta","HTrackEta", 700, -3.5,+3.5); HTrackEta.SetStats(true);
	TH1D HTrackPhi("HTrackPhi","HTrackPhi", 1000, -Pi,+Pi); HTrackPhi.SetStats(true);
	TH1D HTrackInvPt("HTrackInvPt","HTrackInvPt", 1000, -1/2.,+1/2.); HTrackInvPt.SetStats(true);
	TH2I HTrackInvPtvsPhi("HTrackInvPtvsPhi","HTrackInvPtvsPhi",1000, -Pi,+Pi,1000, -1/2.,+1/2.);HTrackInvPtvsPhi.SetStats(true);
	
	TH1D HTrackPz("HTrackPz","HTrackPz", 10000, -100., +100.); HTrackPz.SetStats(true);
	TH1D HTrackInvPz("HTrackInvPz","HTrackInvPz", 10000, -5., +5.); HTrackInvPz.SetStats(true);
	TH1D HTrackPhi0("HTrackPhi0","HTrackPhi0", 1000, -Pi,+Pi); HTrackPhi0.SetStats(true);
	TH2I HTrackInvPzvsPhi0("HTrackInvPzvsPhi0","HTrackInvPzvsPhi0",1000, -Pi,+Pi,1000, -1/2.,+1/2.);HTrackInvPzvsPhi0.SetStats(true);
	
	
	
	
	// Define vector of 1-Dim histograms for barrel hits
	
	const static unsigned nBarrels = dg.nBarrels;
	
	std::vector<TH1D*> HitBX;	
	for(unsigned iB = 0; iB != nBarrels; ++iB){
		stringstream ss;
		ss << "HitBX" << iB;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, 2* dg.B[iB].r * Pi, -dg.B[iB].r * Pi, dg.B[iB].r * Pi); 
		HitBX.push_back(h);	
	}
	
	std::vector<TH1D*> HitBZ;	
	for(unsigned iB = 0; iB != nBarrels; ++iB){
		stringstream ss;
		ss << "HitBZ" << iB;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, -dg.B[iB].zMin + dg.B[iB].zMax, dg.B[iB].zMin, dg.B[iB].zMax); 
		HitBZ.push_back(h);	
	}
	
	std::vector<TH1D*> HitBT;	
	for(unsigned iB = 0; iB != nBarrels; ++iB){
		stringstream ss;
		ss << "HitBT" << iB;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, 2000, 0., 2000.); 
		HitBT.push_back(h);	
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
		
	std::vector<TH1D*> HitDPhi;
	
	for(unsigned iD = 0; iD != nDiscs; ++iD){
		stringstream ss;
		ss << "HitDPhi" << iD;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, 1000, -1., +1.);
		HitDPhi.push_back(h);	
	}
	
	std::vector<TH1D*> HitDR;
	
	for(unsigned iD = 0; iD != nDiscs; ++iD){
		stringstream ss;
		ss << "HitDR" << iD;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, -dg.D[iD].rMin + dg.D[iD].rMax, dg.D[iD].rMin, dg.D[iD].rMax);
		HitDR.push_back(h);	
	}
	
		
	std::vector<TH1D*> HitDT;	
	for(unsigned iD = 0; iD != nDiscs; ++iD){
		stringstream ss;
		ss << "HitDT" << iD;
		string sss = ss.str();
		TString title = TString(sss.c_str());
		TH1D* h = new TH1D(title,title, 4000, 0., 4000.); 
		HitDT.push_back(h);	
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
	
	//const static unsigned nBarrels = dg.nBarrels;
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
	
	//const static unsigned nDiscs = dg.nDiscs;
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
	
	
	// MAIN LOOP ON EVENT GENERATION FOR TRAINING ///////////////////////////////////////////////
	
	
	cout << "BEGIN TRAINING PHASE " << TrainingPhase << endl;
	
	TrackGeometry tgoriginal = tg; // make copy of original track geometry
	
	double phiMin = par.HTA_t_phi - par.HTA_t_deltaPhi;
	double etaMin = par.HTA_t_eta - par.HTA_t_deltaEta;
	double invPtmin = par.HTA_t_invPt_min;
	
	tg.t_deltaPhi /= (double)HTA_NphiBins;
	tg.t_deltaEta /= (double)HTA_NetaBins;
	
	phiMin += tg.t_deltaPhi;
	etaMin += tg.t_deltaEta;
	
	double deltaInvPt = (par.HTA_t_invPt_max - par.HTA_t_invPt_min)/(double)HTA_NinvptBins;
	
	double minPhi = 0., maxPhi = 0.;
	
// LOOP ON ALL CELLS OF HTM ARRAY ////////////////////////////////////////

	 for(int i = 0; i != HTA_NphiBins; ++i){ 
	 		cerr << " " << i;
			for(int j = 0; j != HTA_NetaBins; ++j)
				for(int k = 0; k != HTA_NinvptBins; ++k){			

					//cerr << "Training cell [" <<i<<"]["<<j<<"]["<<k<<"]"<< endl;
			
		
					///// modify tg to generate tracks only inside cell
		
					tg.t_phi = phiMin + tg.t_deltaPhi*2.*(double)i;
					tg.t_eta = etaMin + tg.t_deltaEta*2.*(double)j;
					tg.t_invPt_min = par.HTA_t_invPt_min + (double)k*deltaInvPt;
					tg.t_invPt_max = tg.t_invPt_min + deltaInvPt;

			
					Event ev(dg, nTracks);// one event with nTracks tracks	
	
					// LOOP ON TRACKS
		
					unsigned nTracks = ev.trackList.size();
					for(unsigned iT = 0; iT != nTracks; ++iT) { // begin loop on tracks
		
						Track thisTrack = ev.trackList[iT];	
			
				
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
	
						// loop on hits of each track
								
						unsigned nHitB = 0;
						unsigned nHitD = 0;	
						unsigned nHits = thisTrack.hitList.size();
		
						for(unsigned iH = 0; iH != nHits; ++iH) { // begin loop on hits of each track

							Hit thisHit = thisTrack.hitList[iH]; /// This hit
				
							// Train HT Array	/////////////////////////////////////////			
				
							HTA.train(thisHit,thisTrack,i,j,k);
				
							/////////////////////////////////////////////////////////////
										
							double tx = thisHit.timeExpected(dg, Mass[2], thisTrack.invPt);
				
							if(thisHit.hitType == 'B') {
								++nHitB;	
								HitBTX[thisHit.iLayer]->Fill(thisHit.t - tx);
							}
										
							if(thisHit.hitType == 'D') {
								++nHitD;
								HitDTX[thisHit.iLayer]->Fill(thisHit.t - tx);
							}
						
						} // end loop on hits
			
						double barrelFraction = (double)nHitB/((double)nHitB + (double)nHitD);
						HBarrelFraction.Fill(barrelFraction);
						HNDvsNB.Fill(nHitB,nHitD);

			
					} // end loop on tracks
		
		
					// LOOP ON HITS
	
					unsigned nHits = ev.hitList.size();
					for(unsigned iH = 0; iH != nHits; ++iH) {
		
						Hit thisHit = ev.hitList[iH];
						double X, Y, Z;
						thisHit.XYZ(dg, X, Y, Z);
						double R = sqrt(X*X + Y*Y);
						double Phi = atan2(Y,X);
			
						if(Phi < minPhi) minPhi = Phi;
						if(Phi > maxPhi) maxPhi = Phi;
						
					/*	if(Phi < - Pi/2.) {
							cout << "******* Phi is " << Phi << endl;
							thisHit.print(cout);
							cout << "X = " << X << " Y = " << Y << " Z = " << Z << endl;
							cout << "R = " << R << " Phi = " << Phi << endl;
							
							//return 0;
						}*/
			
						HitX.Fill(X);
						HitY.Fill(Y);
						HitZ.Fill(Z);
			
						HitZPhi.Fill(Z,Phi);
						HitRPhi.Fill(R,Phi);
			
						if(thisHit.hitType == 'B') {
							HitXYbarrel.Fill(X,Y);
							HitRZ.Fill(Z,R);
							HitBXZ[thisHit.iLayer]->Fill(thisHit.x1,Z);
							HitBX[thisHit.iLayer]->Fill(thisHit.x1);
							HitBZ[thisHit.iLayer]->Fill(thisHit.x2);
							HitBT[thisHit.iLayer]->Fill(thisHit.t);
						}
						if(thisHit.hitType == 'D') {
							HitXYdisc.Fill(X,Y);
							HitRZ.Fill(Z,R);
							HitDPhiR[thisHit.iLayer]->Fill(Phi,R);
							HitDPhi[thisHit.iLayer]->Fill(Phi);
							HitDR[thisHit.iLayer]->Fill(R);
							HitDT[thisHit.iLayer]->Fill(thisHit.t);
						}	
					}// end loop on hits
	
				}
	} // end loop on HTM array cells ///////////////////////////////////////////////////////////
	
		
	cout << "Event generation complete" << endl;

	
	cout << endl;
	cout << "Phi limits for Bib distribution" << endl;
	cout << "(copy and paste into Parameters.h)"<< endl;
	cout << " double gen_phia = " << minPhi << ";" << endl;
	cout << " double gen_phib = " << maxPhi << ";" << endl;
	cout << endl;
		
	
	if(TrainingPhase==1 && Diagonalize) {
		cout << "Diagonalize ...";
		HTA.diagonalize();
		cout << endl;
	}
		
		
	
	if(TrainingPhase == 1 || TrainingPhase == 2) {
		// write data file
		cout << "Writing data file..." << endl;
		ofstream outfile;	
		outfile.open(dataFileName);
		HTA.write(outfile);
		outfile.close();
	}
	
	if(TrainingPhase != 3 ) {				
		cout << "Write histogram file " << train_histFileName.c_str() << " ...";
		histFile->Write(); // write histogram file
		cout << endl;
		return 0;// end main without training results
	}
	
	// This is executed only for TrainingPhase == 3 
	
	
	//if(Summary) HTA.print(cout,1);
	
	unsigned int nPhi = HTA.NphiBins;
	unsigned int nEta = HTA.NetaBins;
	unsigned int nInvpt = HTA.NinvptBins;
	//unsigned int nLayers = dg.nBarrels + dg.nDiscs;
	
	// Loop on all cells of HT array (3 nested loops)
	for(unsigned int iPhi = 0; iPhi != nPhi; ++iPhi)
		for(unsigned int iEta = 0; iEta != nEta; ++iEta)
			for(unsigned int iInvpt = 0; iInvpt != nInvpt; ++iInvpt){		
				HTArrayElement thisElement = HTA.ArrElem[iPhi][iEta][iInvpt];
				unsigned int kLayers = thisElement.layerIndHitStat.size();			
				for(auto it = thisElement.layerIndHitStat.begin(); it != thisElement.layerIndHitStat.end(); ++it){
					unsigned int N = it->second.nEntries;
        			int layerInd = it->first;
        			double d0 = it->second.u0Max - it->second.u0Min;
        			double d1 = it->second.u1Max - it->second.u1Min;
        			double d2 = it->second.u2Max - it->second.u2Min;
        			Hd0Width.Fill(d0); 
        			Hd1Width.Fill(d1);
        			Hd2Width.Fill(d2);
        		}	
			
			} // end loop on cells of HT array for summarizing
			
			
	if(!Summary) {			
		cout << "Write histogram file " << train_histFileName.c_str() << " ...";
		histFile->Write(); // write histogram file
		cout << endl;
		return 0;// end main without summary
	}
	
	//this is executed only if TrainingPhase == 3 && Summary == true		

	
	// summarize coordinate cut statistics
	cout << "Summarizing training results" << endl;
		
	// Loop on all cells of HT array (3 nested loops)
	for(unsigned int iPhi = 0; iPhi != nPhi; ++iPhi)
		for(unsigned int iEta = 0; iEta != nEta; ++iEta)
			for(unsigned int iInvpt = 0; iInvpt != nInvpt; ++iInvpt){		
				HTArrayElement thisElement = HTA.ArrElem[iPhi][iEta][iInvpt];
				unsigned int kLayers = thisElement.layerIndHitStat.size();
				cout << endl;
				cout << "iPhi: " << iPhi << " iEta: " << iEta << " iInvpt: " << iInvpt << " kLayers: " << kLayers << endl ;		
	 
				for(auto it = thisElement.layerIndHitStat.begin(); it != thisElement.layerIndHitStat.end(); ++it){
					unsigned int N = it->second.nEntries;
        			int layerInd = it->first;
        			double d0 = it->second.u0Max - it->second.u0Min;
        			double d1 = it->second.u1Max - it->second.u1Min;
        			double d2 = it->second.u2Max - it->second.u2Min;
        		    cout << "layer: " << layerInd;			
					cout << " N: " << N;
					cout << " d0: " << d0;
					cout << " d1: " << d1;
					cout << " d2: " << d2;
					cout << endl;
        		}	
			
			} // end loop on cells of HT array for summarizing
			
		
	cout << "Write histogram file " << train_histFileName.c_str() << " ...";
	histFile->Write(); // write histogram file
	cout << endl;

			
	return 0;		
	
} // end main



