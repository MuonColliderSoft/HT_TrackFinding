//
//  Hit.h
//  MuonColliderToy
//
//  Created by Luciano Ristori on 2/28/23
//


    class Hit {
    
    public:
    
        unsigned ID; // sequence number in hitList contained in event
        char hitType; // 'B' for barrel, 'D' for disc, 'G' for ghost (other codes possible)       
        unsigned iLayer; // layer number (separate count of barrels and layers)
        unsigned layerInd; // single index going through all barrels and discs
        double x1; // primary coordinate mm
        double x2;  // secondary coordinate mm
        double t; //  time coordinate  mm 
        double u0; 
        double u1;
        double u2; 
        unsigned trackInd; // index of parent track in track list (member of Event)
               
        Hit(char _hitType, int _iLayer, double _x1, double _x2, double t, unsigned trackInd);
        Hit(char _hitType); 
        void XYZ(DetectorGeometry &g, double &X, double &Y, double &Z);  
		void print(ostream &out); 
		void write(ostream &out);
		void printXYZ(DetectorGeometry &g, ostream &out);
		double timeExpected(DetectorGeometry &g, double mass, double invPt);
		double realPhi(DetectorGeometry &g);
		// TOF of a particle of mass "mass" and inverse pt "invPt" from the origin to the hit point
		bool isSeed();
    
    };
    
   
        
      
    