//
//  BibFileReader.cpp
//  HT_TrackFinding
//
//  Created by Luciano Ristori on 1/31/23
//

 //////////////////////////////////////////////////////    

    
   class Record {
   
	public:
	 	char hitType ; // B (barrel) or D (disc)
 	 	unsigned iLayer ; // indexed separately for barrels and discs
 	 	double x2 ; // secondary coordinate (mm)
 	 	double time; //  time (mm)
 	//------------------------------------------------- 
	void print(std::ostream &out){
		out << "type: " << hitType << " iLayer: " << iLayer << " x2: " << x2 << " time: " << time << std::endl;
	} 
	//-------------------------------------------------
	void write(std::ostream &out){
		out << hitType << " " << iLayer << " " << x2 << " " << time << std::endl;
	} 
	//-------------------------------------------------
	void read(std::istream &in){
		in >> hitType >> iLayer >> x2 >> time;
	} 
};

 //////////////////////////////////////////////////////    
    
    class BibFileReader {
        
    public:
        
        std::vector<Record> bibList;
        unsigned nBIB = 0;
        
        double phiMin;
        double phiMax;
        
        
		//-------------------------------------------------
        
        void setPhiLimits(double _phiMin, double _phiMax){
        	phiMin = _phiMin;
        	phiMax = _phiMax;
        }
        
        //-------------------------------------------------
        int readFile(std::string BibFileName) {
        
       		// Open data file
	
			std::cout << "Opening Data File: " << BibFileName << " ..." << std::endl;
			std::ifstream infile;
			infile.open(BibFileName);
			Record dummyRecord;
	
			if(infile){ // data file found
				for(; infile ;){
					dummyRecord.read(infile);
					bibList.push_back(dummyRecord);
				}
				nBIB = bibList.size();
				return 0;		
			} 
	
			else {
				std::cout << "Data File not found, abort." << std::endl;
				return -1;
			}
	     
        }
        
        //-------------------------------------------------
        unsigned size(){
        	return nBIB;
        }
        
        //-------------------------------------------------      
        Hit randomBibHit(){
        	std::uniform_int_distribution<unsigned> intDistribution(0,nBIB-1);
        	unsigned bibIndex = intDistribution(generator);
        	Record rec = bibList[bibIndex]; 
        	double x1Min = 0.; double x1Max = 0.;
			if(rec.hitType == 'B') {
				x1Max = phiMax*g.B[rec.iLayer].r;
				x1Min = phiMin*g.B[rec.iLayer].r;
        	}
        	else{
				x1Max = phiMax*rec.x2;
				x1Min = phiMin*rec.x2;
        	}
        	double alfa = distribution(generator); 	
        	double x1 =  x1Min + (x1Max - x1Min) * alfa; 
        	//cout << x1Min << " " << x1Max << " " <<  alfa << endl;      
        	Hit bibHit(rec.hitType, rec.iLayer , x1, rec.x2, rec.time, 0);        	
        	return bibHit;
        
        }
        
        
    }; // end class BibFileReader
    
