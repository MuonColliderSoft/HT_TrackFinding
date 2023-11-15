//
//  DetectorGeometry.cpp
//  MuonColliderToy
//
//  Created by Luciano Ristori on 64/28/23.
//


extern Parameters par;


class Barrel {

//
// Ideal cylindrical detector along z

public:

	// detector boundaries
	
	double r;
	double zMin;
	double zMax;
	
	// gaussian measurement precisions (sigma)
		
	double zPrec; // along z
	double xphiPrec; // linear coordinate in phi direction
	double tPrec; // time 
	
	
	void print(std::ostream &out) {
	
		out << "r = " << r << " zMin = " << zMin << " zMax = " << zMax 
			<< " xphiPrec = " << xphiPrec << " zPrec = " << zPrec << " tPrec = " << tPrec << endl;	
			
	}// end print

}; // end Barrel




class Disc {

public:

	// detector boundaries

	double z;
	double rMin;
	double rMax;
		
	// gaussian measurement precisions (sigma)
	
	double rPrec; // along r
	double xphiPrec; // linear coordinate in phi direction
	double tPrec; // time
	
	
	
	void print(std::ostream &out) {
	
		out << "z = " << z << " rMin = " << rMin << " rMax = " << rMax 
			<< " xphiPrec = " << xphiPrec << " rPrec = " << rPrec << " tPrec = " << tPrec << endl;	
	
	}// end print
	
	
}; // end Disc


class DetectorGeometry {


public:

	double xphiSigmaB;
	double xphiSigmaD;
	double zSigmaB;
	double rSigmaD; 
	double tSigmaB;
	double tSigmaD;	
   

	// DETECTORS	

	const static unsigned nBarrels = 11;
	const static unsigned nDiscs = 30;	

	Barrel _B, B[nBarrels];
	Disc _D, D[nDiscs];
	
	// one barrel used only as boudary for all tracks
	const static unsigned iBoundaryBarrel = nBarrels - 1;

	
	DetectorGeometry() { // default constructor
	
		// detector measurement errors
		
		xphiSigmaB = par.geo_xphiSigmaB;
		xphiSigmaD = par.geo_xphiSigmaD;
		zSigmaB = par.geo_zSigmaB;
		rSigmaD = par.geo_rSigmaD; 
		tSigmaB = par.geo_tSigmaB;
		tSigmaD = par.geo_tSigmaD;
			

		unsigned iB = 0;
		_B.r = 31.; _B.zMax = 65.; _B.zMin =-_B.zMax; _B.xphiPrec = xphiSigmaB; _B.zPrec = zSigmaB; _B.tPrec = tSigmaB;    B[iB++] = _B;
		_B.r = 52.; _B.zMax = 65.; _B.zMin =-_B.zMax; _B.xphiPrec = xphiSigmaB; _B.zPrec = zSigmaB; _B.tPrec = tSigmaB;    B[iB++] = _B;
		_B.r = 75.; _B.zMax = 65.; _B.zMin =-_B.zMax; _B.xphiPrec = xphiSigmaB; _B.zPrec = zSigmaB; _B.tPrec = tSigmaB;    B[iB++] = _B;
		_B.r = 103.; _B.zMax = 65.; _B.zMin =-_B.zMax; _B.xphiPrec = xphiSigmaB; _B.zPrec = zSigmaB; _B.tPrec = tSigmaB;    B[iB++] = _B;
		_B.r = 129.; _B.zMax = 480.; _B.zMin =-_B.zMax; _B.xphiPrec = xphiSigmaB; _B.zPrec = zSigmaB; _B.tPrec = tSigmaB;    B[iB++] = _B;
		_B.r = 350.; _B.zMax = 480.; _B.zMin =-_B.zMax; _B.xphiPrec = xphiSigmaB; _B.zPrec = zSigmaB; _B.tPrec = tSigmaB;    B[iB++] = _B;
		_B.r = 550.; _B.zMax = 700.; _B.zMin =-_B.zMax; _B.xphiPrec = xphiSigmaB; _B.zPrec = zSigmaB; _B.tPrec = tSigmaB;    B[iB++] = _B;
		_B.r = 810; _B.zMax = 1260.; _B.zMin =-_B.zMax; _B.xphiPrec = xphiSigmaB; _B.zPrec = zSigmaB; _B.tPrec = tSigmaB;    B[iB++] = _B;
		_B.r = 1150.; _B.zMax = 1260.; _B.zMin =-_B.zMax; _B.xphiPrec = xphiSigmaB; _B.zPrec = zSigmaB; _B.tPrec = tSigmaB;    B[iB++] = _B;
		_B.r = 1510.; _B.zMax = 1260.; _B.zMin =-_B.zMax; _B.xphiPrec = xphiSigmaB; _B.zPrec = zSigmaB; _B.tPrec = tSigmaB;    B[iB++] = _B;
		
		_B.r = 1511.; _B.zMax = 5000.; _B.zMin =-_B.zMax; _B.xphiPrec = xphiSigmaB; _B.zPrec = zSigmaB; _B.tPrec = tSigmaB;    B[iB++] = _B;// Boundary
		
		
		
		unsigned iD = 0;
		_D.z = 80.; _D.rMin = 27.; _D.rMax = 113.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = 120.; _D.rMin = 31.; _D.rMax = 113.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = 200.; _D.rMin = 40.; _D.rMax = 113.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = 280.; _D.rMin = 52.; _D.rMax = 113.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = 520.; _D.rMin = 100.; _D.rMax = 440.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = 800.; _D.rMin = 150.; _D.rMax = 560.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = 1110.; _D.rMin = 200.; _D.rMax = 560.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = 1380.; _D.rMin = 230.; _D.rMax = 560.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = 1660.; _D.rMin = 250.; _D.rMax = 560.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = 1950.; _D.rMin = 270.; _D.rMax = 560.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = 2200.; _D.rMin = 280; _D.rMax = 560.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = 1300.; _D.rMin = 620.; _D.rMax = 1440.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = 1610.; _D.rMin = 620.; _D.rMax = 1440.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = 1890.; _D.rMin = 620.; _D.rMax = 1440.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = 2200.; _D.rMin = 620.; _D.rMax = 1440.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = -80.; _D.rMin = 27.; _D.rMax = 113.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = -120.; _D.rMin = 31.; _D.rMax = 113.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = -200.; _D.rMin = 40.; _D.rMax = 113.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = -280.; _D.rMin = 52.; _D.rMax = 113.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = -520.; _D.rMin = 100.; _D.rMax = 440.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = -800.; _D.rMin = 150.; _D.rMax = 560.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = -1110.; _D.rMin = 200.; _D.rMax = 560.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = -1380.; _D.rMin = 230.; _D.rMax = 560.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = -1660.; _D.rMin = 250.; _D.rMax = 560.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = -1950.; _D.rMin = 270.; _D.rMax = 560.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = -2200.; _D.rMin = 280; _D.rMax = 560.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = -1300.; _D.rMin = 620.; _D.rMax = 1440.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = -1610.; _D.rMin = 620.; _D.rMax = 1440.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = -1890.; _D.rMin = 620.; _D.rMax = 1440.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;
		_D.z = -2200.; _D.rMin = 620.; _D.rMax = 1440.; _D.xphiPrec = xphiSigmaD; _D.rPrec = rSigmaD; _D.tPrec = tSigmaD;    D[iD++] = _D;


		cout << " DETECTOR GEOMETRY INITIALIZED" << endl;
		if(iB == nBarrels && iD == nDiscs) return;
		else cout << " DETECTOR GEOMETRY INITIALIZATION ERROR" << endl;
		return;
		
		
	}; // end default constructor
	
	
	
	
	void print(std::ostream &out) {
	
		out << endl << "DETECTOR GEOMETRY DATA" << endl;
		
		out << endl;
		out  << "Barrels: "<< nBarrels << endl;
		for(unsigned ib = 0; ib != nBarrels; ++ ib){
			out << "L" << ib << "   Barrel " << ib  << ": " ;
			B[ib].print(out);	
		}
		out << endl;
		
		out  << "Discs: " << nDiscs << endl;
		for(unsigned id = 0; id != nDiscs; ++ id){
			out << "L" << id + nBarrels << "  Disc " << id  << ": ";
			D[id].print(out);	
		}
		
		out << endl;
		
	}// end print


}; // end class DetectorGeometry





