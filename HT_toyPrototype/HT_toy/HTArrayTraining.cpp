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

std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.,1.);
std::normal_distribution<double> gauss(0.0,1.0);

#include "GlobalConstants.cpp"

#include "Parameters.h"
Parameters par; // class that holds all configurable parameters					   

//#include "Geometry.cpp"
//Geometry g;

#include "DetectorGeometry.cpp"
DetectorGeometry dg;

#include "TrackGeometry.cpp"
TrackGeometry tg;

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
					   
unsigned nEvents; // number of events to be generated for training		   
bool Diagonalize;
bool Summary;
bool Special;



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



int main(){
  

	// Parameter initialization:

	nEvents = par.train_nEvents;
	Diagonalize = par.train_Diagonalize;
	Summary = par.train_Summary;
	

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
	
	HTA.initHist3D();
  

  
	dg.print(cout);
	
	
	TH1D HBarrelFraction("HBarrelFraction","HBarrelFraction", 11, -0.05, 1.05);HBarrelFraction.SetStats(false);
	TH2D HNDvsNB("NDvsNB","NDvsNB",  13, -0.5,+12.5, 13, -0.5, +12.5); HNDvsNB.SetStats(false);
	
	
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
	

    for(unsigned iEv = 0; iEv != nEvents; ++iEv){
    
    	//cout << "New Event " << iEv << endl;
    	
    	
    	if(nEvents >= 1e6)
    		if(iEv && ((iEv % (int)1e5) == 0) ) cout << iEv << "/" << nEvents << " processed events" << endl;
    
    		
		Event ev(dg, 1);// one event with 1 track
		
		//ev.print(cout,1);
		
	
	
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
				
				HTA.train(thisHit,thisTrack);
				
				/////////////////////////////////////////////////////////////
										
				double tx = thisHit.timeExpected(dg, mass[2], thisTrack.invPt);
				
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
			
			
			HitX.Fill(X);
			HitY.Fill(Y);
			HitZ.Fill(Z);
			
			//if(Z < 0.) cout << "iEv = " << iEv << " iH = " << iH << " layer = " << thisHit.iLayer << " Z = " << Z << endl;
			
			HitZPhi.Fill(Z,Phi);
			HitRPhi.Fill(R,Phi);
			
			if(thisHit.hitType == 'B') {
				HitXYbarrel.Fill(X,Y);
				HitRZ.Fill(Z,R);
				HitBXZ[thisHit.iLayer]->Fill(thisHit.x1,Z);
				HitBX[thisHit.iLayer]->Fill(thisHit.x1);
				HitBZ[thisHit.iLayer]->Fill(thisHit.x2);
				HitBT[thisHit.iLayer]->Fill(thisHit.t);
				//double tx = thisHit.timeExpected(g,mass[1],)/////////////////////////
				
				//double xCorr = thisHit.x1 - R*(dg.phi0Center + 6.e-4 * dg.invPzCenter * Z);
				//HitBXCorr[thisHit.iLayer]->Fill(xCorr);
			}
			if(thisHit.hitType == 'D') {
				HitXYdisc.Fill(X,Y);
				HitRZ.Fill(Z,R);
				HitDPhiR[thisHit.iLayer]->Fill(Phi,R);
				HitDPhi[thisHit.iLayer]->Fill(Phi);
				HitDR[thisHit.iLayer]->Fill(R);
				HitDT[thisHit.iLayer]->Fill(thisHit.t);
					//if(thisHit.iLayer==10)cout << thisHit.t << "  ";////////////////////
				
				//double phiCorr = Phi - dg.t_phi - asin(6.e-4*dg.t_invPt_mean*R);
				//HitDPhiCorr[thisHit.iLayer]->Fill(phiCorr);
			}	
		}// end loop on hits
		
		//ev.print(cout,1); // mode = 1 prints all hits for each track
		//ev.printXYZ(g,cout); // prints data for Mathematica Plot
	
	} // end loop on events ///////////////////////////////////////////////////////////
	
	
	cout << "Event generation complete" << endl;
	
	
	cout << " Write histogram file " << train_histFileName.c_str() << " ...";
	histFile->Write(); // write histogram file
	cout << endl;
	
	
	
	if(TrainingPhase==1 && Diagonalize) {
		cout << "Diagonalize ...";
		HTA.diagonalize();
		cout << endl;
	}
		if(Summary) HTA.print(cout,1);
		
	
	if(TrainingPhase == 1 || TrainingPhase == 2) {
		// write data file
		cout << "Writing data file..." << endl;
		ofstream outfile;	
		outfile.open(dataFileName);
		HTA.write(outfile);
		outfile.close();
	}
	
	
	
	
	if(!Summary) return 0; // end main without training results
	if(TrainingPhase != 3 ) return 0;// end main without training results
	
	
	
	
	// summarize coordinate cut statistics
	cout << "Summarizing training results" << endl;
	
	unsigned int nPhi = HTA.NphiBins;
	unsigned int nEta = HTA.NetaBins;
	unsigned int nInvpt = HTA.NinvptBins;
	unsigned int nLayers = dg.nBarrels + dg.nDiscs;
	
	
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
					//if(N < 100) continue; // do not consider cell layers with low statistics
        			int layerInd = it->first;
        			double d0 = it->second.u0Max - it->second.u0Min;
        			double d1 = it->second.u1Max - it->second.u1Min;
        			double d2 = it->second.u2Max - it->second.u2Min;
        			//double d0 = sqrt(it->second.u0v);
        			//double d1 = sqrt(it->second.u1v);
        			//double d2 = sqrt(it->second.u2v);
        		    cout << "layer: " << layerInd;			
					cout << " N: " << N;
					cout << " d0: " << d0;
					cout << " d1: " << d1;
					cout << " d2: " << d2;
					cout << endl;
        		}	
			
			} // end loop on cells of HT array for summarizing
			
			
			
			
	return 0;		
	
} // end main



