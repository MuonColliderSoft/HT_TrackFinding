
 #ifndef MDHT_BIBFILEREADER_CPP
 #define MDHT_BIBFILEREADER_CPP
 
//////////////////////////////////////////////////////    
//
//  BibFileReader.cpp
//  HT_TrackFinding
//
//  Created by Luciano Ristori on 1/31/23
//
//  Modified for NewHTATraining - by Luciano 9/25/25
//
//////////////////////////////////////////////////////    
    
    class BibFileReader {
        
    public:
    
    	std::vector<Hit> bibList;
    
    	
    	
    	int readFile(std::string BibFileName){
    	    			
			TrackReader* bibReader = new TrackReader(BibFileName);
			Track thisTrack;
			
			bool eof = bibReader->read(thisTrack);
        		if(eof) {
        			cerr << "BibFileReader failed" << endl;
        			
        		}		
        		// loop on hits of this track
        		     		     	
        		for(int iH = 0; iH != (int)thisTrack.hitList.size(); ++iH){
        			Hit thisHit = thisTrack.hitList[iH];
        			thisHit.trackInd = 0; // background hit	
        			bibList.push_back(thisHit); // add this hit to bibList
        			
        		} // end loop on hits of this track
        		        		
        		cout << "Created BIB Pool " << endl;
        		return 0;
        		
    	}
    	
        
        //-------------------------------------------------      
         Hit randomBibHit(){
        	std::uniform_int_distribution<unsigned> intDistribution(0,bibList.size()-1);
        	unsigned bibIndex = intDistribution(generator);
        	return bibList[bibIndex];
        }
        
        int size(){
        	return (int)bibList.size();
        }
        
        
    }; // end class BibFileReader
    
#endif // MDHT_BIBFILEREADER_CPP
    
