//
//  Event.cpp
//  MuonColliderToy
//
//  Created by Luciano Ristori on 6/30/21
//


    
   
    
    
    
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
        			double hitPhi, hitZ, hitT;
        			if(thisTrack.phizBarrel(g, iB, hitPhi, hitZ, hitT)){ // find hit	
        			 Hit thisHit('B', iB, hitPhi, hitZ, hitT, iT);
        			 if(!thisHit.isSeed())continue;
        			 	thisHit.ID = hitList.size();
        			 	hitList.push_back(thisHit);
        			 	thisTrack.hitList.push_back(thisHit); // add copy of this hit in generating track     			
        			}	
        		}// end loop ond barrels
        		
        		// find hits in all disc detectors
        		for(unsigned iD = 0; iD != g.nDiscs; ++iD){
        			double hitPhi, hitR, hitT;
        			if(thisTrack.xyDisc(g, iD, hitPhi, hitR, hitT)){ // find hit		
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
        
        void addBibHits(BibFileReader &bibRead, unsigned nBibHits){
        	for(unsigned iH = 0; iH != nBibHits; ++iH){
        		Hit bibHit = bibRead.randomBibHit();
        		bibHit.ID = hitList.size();
        		hitList.push_back(bibHit);
        	}
        }
        
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
    
