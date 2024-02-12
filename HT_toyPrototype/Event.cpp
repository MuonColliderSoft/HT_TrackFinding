//
//  Event.cpp
//  MuonColliderToy
//
//  Created by Luciano Ristori on 6/30/21
//


    extern bool verbose;
    
    
    class Event {
        
    public:
        
        vector<Track> trackList;
        vector<Hit> hitList;
        
        // Constructor
        
        Event(DetectorGeometry &g, int nTracks){
        
        	// Create all tracks
        	for(unsigned iT = 1; iT != nTracks+1; ++iT){// iT=0 is BIB
        		
        		Track thisTrack(tg); //  create track
        		thisTrack.ID = iT; // assign track iD = track index
        		
        		// find hits in all barrel detectors
        		for(unsigned iB = 0; iB != g.nBarrels; ++iB){
        			if(iB == g.iBoundaryBarrel) continue; // skip boundary barrel
        			if(distribution(generator_trk) < par.geo_ineffB)continue;
        			double hitPhi, hitZ, hitT;
        			bool smear = true;
        			if(thisTrack.phizBarrel(g, iB, hitPhi, hitZ, hitT, smear)){ // find hit	
        			 Hit thisHit('B', iB, hitPhi, hitZ, hitT, iT);
        			 if(!thisHit.isSeed())continue;
        			 	thisHit.ID = hitList.size();
        			 	hitList.push_back(thisHit);
        			 	thisTrack.hitList.push_back(thisHit); // add copy of this hit in generating track     			
        			}	
        		}// end loop ond barrels
        		
        		// find hits in all disc detectors
        		
        		for(unsigned iD = 0; iD != g.nDiscs; ++iD){
        		
        			double hitPhi, hitZ, hitT;
        			bool smear = false;
        			if(thisTrack.phizBarrel(g, g.iBoundaryBarrel, hitPhi, hitZ, hitT, smear)){ // find hit in boundary barrel
        				if(fabs(hitZ) < fabs(g.D[iD].z)) continue;
        			}
        				
        			if(distribution(generator_trk) < par.geo_ineffD)continue;
        			
        			double hitR;       			
        			smear = true;
        			if(thisTrack.xyDisc(g, iD, hitPhi, hitR, hitT, smear)){ // find hit		
        			 Hit thisHit('D', iD, hitPhi, hitR, hitT, iT);
        			 if(!thisHit.isSeed())continue;
        			 	thisHit.ID = hitList.size();
        			 	hitList.push_back(thisHit);
        			 	thisTrack.hitList.push_back(thisHit); // add copy of this hit in generating track     			
        			}	
        		}// end loop on discs
        		
        		trackList.push_back(thisTrack); // add this track to trackList
        		
        	}// end loop on tracks
        		
        	
        }// end event constructor
        
        
        
        // Adds BIB hits to the event
        
        void addBibHits(BibFileReader &bibRead, unsigned nBibHits){
        
        	for(unsigned iH = 0; iH != nBibHits; ++iH){
        		Hit bibHit = bibRead.randomBibHit();
        		bibHit.ID = hitList.size();
        		hitList.push_back(bibHit);
        	}
        	
        	// shuffle hits
        	
  			std::shuffle(std::begin(hitList), std::end(hitList), generator);
  
  			if(verbose) cout << "hitList.size():" << hitList.size() << endl;
  			
  			// redefine hit ID's
  			
  			unsigned IDmap[100];// max possible size 		
  			
  			for(unsigned iH = 0; iH != hitList.size(); ++iH){
  				
  				Hit thisHit = hitList[iH];
  				unsigned trackInd = thisHit.trackInd;				
  				if(trackInd != 0) { // not BIB hit
  				
  					// create IDmap from old hit ID to new hit ID
  											
  					unsigned IDoldmax = 0;
					unsigned IDold = thisHit.ID;
					//cout << "IDold: " << IDold << endl;
					hitList[iH].ID = iH; // set new ID in hitList
					if(IDold < 100) IDmap[IDold] = iH;
					else cout << "*** ERROR *** IDmap overflow " << endl;
					if(IDold > IDoldmax) IDoldmax = IDold;
					//cout << "IDoldmax: " << IDoldmax << endl;
				}
							
				hitList[iH].ID = iH; // update to new position in list
				
  			} // end loop on hits
  						
  			// set new ID in hitList inside all tracks using IDmap
  			
  			for(unsigned iT = 0; iT != trackList.size(); ++iT){
  				for(unsigned iH = 0; iH != trackList[iT].hitList.size(); ++iH){
  					trackList[iT].hitList[iH].ID = IDmap[trackList[iT].hitList[iH].ID];				
  				}		 			
  			}
  			  			
  			return;
  			
        }// end addBibHits
        
        
        
        void print(ostream &out, int mode = 0){
        
        	// mode = 0: print only track parameters
        	// mode = 1: print all hits
        
        	for(unsigned iT = 0; iT != trackList.size(); ++iT) {
        		trackList[iT].print(out,mode);
        	//	if(mode == 0) continue;
        	//	unsigned nHits = trackList[iT].hitList.size();
        	//	for(unsigned iH = 0; iH != nHits; ++iH) trackList[iT].hitList[iH].print(out); 		
        	}
        	
        }// end print
        
        
         void printXYZ(DetectorGeometry &g, ostream &out){
         
        // Print input for Mathematica ListPointPlot3D   
        
        out << "{";
        	for(unsigned iH = 0; iH != hitList.size(); ++iH) {
        		hitList[iH].printXYZ(g, out);
        		if(iH != hitList.size() - 1) out << ",";
        		else out << "}";
        		out << endl;
        	}
        	
        }// end printXYZ
        
        
        
        
    }; // end class Event
    
