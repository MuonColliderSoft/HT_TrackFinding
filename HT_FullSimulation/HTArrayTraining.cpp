#include <iostream> 
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <random>

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"


using namespace std;

#include "Hit.h"
#include "Hit.cpp"

#include "Track.h"
#include "Track.cpp"
#include "Parameters.h"
#include "TrackReader.h"
#include "HTAmapper.h"

#include "HTArray.cpp"

// random generator

//std::default_random_engine generator;
std::mt19937 generator; // Mersenne Twister 
std::mt19937 generator_trk; // Mersenne Twister for track generation
std::uniform_real_distribution<double> distribution(0.,1.);
std::normal_distribution<double> gauss(0.0,1.0);


Parameters par;

int TrainingPhase;



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



int main(){

	int start = 0;
	
	
	
	vector<Track> trackLists[par.HTA_NphiBins][par.HTA_NetaBins][par.HTA_NinvptBins];

	// Instantiate Hough Transform Array
	
	cout << "Instantiating Hough Transform Array..." << endl;
	
	HTArray *HTA = new HTArray();
			
	//cout << "Done" << endl;
	
	HTA->print(cout);

			
	// Open data file
	
	string dataFileName = par.train_dataFileName;
	cout << "Opening Data File: " << dataFileName << " ..." << endl;
	
	

	ifstream infile;
	infile.open(dataFileName);
	
	if(start != 1)	{
	
		if(infile){ // data file found
			int retcode = HTA->read(infile);
			cout << "HTA.read retcode: " << retcode << endl;
			TrainingPhase = retcode;		
		} 
	
		else {
			cout << "Data File not found, starting training from scratch" << endl;
			TrainingPhase = 1;
		}
	}
	
	else TrainingPhase = 1;
	
		
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

	
	HTA->initHists();
	
	
	TH1D HTrackZ0("HTrackZ0","HTrackZ0", 600, -300,+300); HTrackZ0.SetStats(true);
	TH1D HTrackT0("HTrackT0","HTrackT0", 600, -300,+300); HTrackT0.SetStats(true);
	TH1D HTrackEta("HTrackEta","HTrackEta", 700, -3.5,+3.5); HTrackEta.SetStats(true);
	TH1D HTrackPhi("HTrackPhi","HTrackPhi", 1000, -3.15,+3.15); HTrackPhi.SetStats(true);
	TH1D HTrackInvPt("HTrackInvPt","HTrackInvPt", 1000, -1/2.,+1/2.); HTrackInvPt.SetStats(true);
	
	TH1D HnTracksPerCell("HnTracksPerCell","HnTracksPerCell",1000,0.,1000.);
	
	

  
	// instantiate reader and open file
	
	TrackReader newTrainTrack(par.train_inputTrackFileName);
	
	Track thisTrack;
	//HTAmapper mp;

	long unsigned maxTracks = par.train_maxTracks;
	long unsigned iTrack = 0;
	 
	for(;;){ // infinite loop
	
		if(iTrack%100000 == 0) cout << "iTrack " << iTrack << endl;
		if(iTrack%100000 == 0) cerr << ".";
						
		bool eof = newTrainTrack.read(thisTrack); // read next track from file
		if(eof) break; // End Of File
		if(iTrack == maxTracks && iTrack != 0) break; // max number of tracks reached
		
		
		HTrackZ0.Fill(thisTrack.z0);
		HTrackT0.Fill(thisTrack.t0);
		HTrackEta.Fill(thisTrack.eta);
		HTrackPhi.Fill(thisTrack.phi);
		HTrackInvPt.Fill(thisTrack.invPt);

		
		
		int i, j, k;
		int retcode = trackToCell(i, j, k, thisTrack);
		//cout << endl;
		//cout << retcode << " " << i << " " << j << " " << k;
		if(retcode != 0){
		cout << "******** track " << iTrack << " out of range. Skipped ********** " << endl;
		continue;
		};
		//cout << endl;
		
		if(TrainingPhase == 3 && par.train_WriteFiles)trackLists[i][j][k].push_back(thisTrack);
		
		
		if(par.train_printTracks) thisTrack.print(cout,1); // print track parameters and all hits	
	
		HTA->trainStat(i,j,k); // accumulate training statistics
		
			unsigned nHits = thisTrack.hitList.size();
		
						for(unsigned iH = 0; iH != nHits; ++iH) { // begin loop on hits of each track

							Hit thisHit = thisTrack.hitList[iH]; 	
				
							HTA->train(thisHit,thisTrack,i,j,k);
				
			}	
		
		++iTrack;		
	}
	
	cout << iTrack << " tracks processed" << endl;
	
		
	if(TrainingPhase==1 && par.train_Diagonalize) {
		cout << "Diagonalize ...";
		HTA->diagonalize();
		cout << endl;
	}
	
	
	// write data file
	
		cout << "Writing data file..." << endl;
		ofstream outfile;	
		outfile.open(dataFileName);
		HTA->write(outfile);
		outfile.close();
		
	// Training summary
	
		if(par.train_Summary) HTA->print(cout,1);
		
	// write histogram file	
			
		cout << "Write histogram file " << train_histFileName.c_str() << " ...";
		histFile->Write(); // write histogram file
		cout << endl;

		cout << "Completed training phase number " << TrainingPhase 
			<< " data file: " << dataFileName  
			<< " histogram file " << train_histFileName << endl;
		
	
	// write tracks in each cell
	
		if(TrainingPhase == 3 && par.train_WriteFiles){

			cout << "WRITING TRACK FILES" << endl;

			string testFileName = HTA->makeFileName("test",1,2,3);
			cout << testFileName << endl;

			for(int i = 0; i != par.HTA_NphiBins; ++i)
				for(int j = 0; j != par.HTA_NetaBins; ++j)
					for(int k = 0; k != par.HTA_NinvptBins; ++k){

				ofstream outfile;
				string name = HTA->makeFileName("TracksPerCell/",i,j,k);
				cout << i << " " << j << " " << k << " " << name << endl;
				outfile.open(name);
				if(!outfile) cout << "opening error: " << name << endl;
				for(int iTrack = 0; iTrack != (int)trackLists[i][j][k].size(); ++iTrack){
					vector<Track> thisTrackList = trackLists[i][j][k];
					thisTrack = thisTrackList[iTrack];
					thisTrack.ID = iTrack;
					thisTrack.write(outfile);
				}
				outfile.close();
				
				HnTracksPerCell.Fill(trackLists[i][j][k].size());

			}

		} // end if(TrainingPhase == 3)
		
		delete HTA;

	return 0;
	
				
} // end main



