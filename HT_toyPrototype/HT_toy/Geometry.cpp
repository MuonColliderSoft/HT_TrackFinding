//
//  Geometry.cpp
//  MuonColliderToy
//
//  Created by Luciano Ristori on 6/24/21.
//



// GLOBAL VARIABLES

const static double Pi = 3.14159265;
const static double SpeedOfLight = 0.000299792458;// m/ps

	
// PARTICLE MASSES in GeV/c2

static const unsigned nMasses = 5;
double mass[] = {
	0.,
	0.000511, // electron: 1
	0.139570, // charged pion: 2
	0.493677, // charged kaon: 3
	0.938272  // proton: 4
	};

bool verbose = false;






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


class Geometry {


	public:
			
    
    // TRACK GENERATION PARAMETERS
    
	//const double t_eta = 0;
	//const double t_deltaEta = 3.;
	
	//const double t_phi = 0.;
	//const double t_deltaPhi = Pi;
	
	//const double t_invPt_max = 1./3.; // Gev/c^(-1)  
	//const double t_invPt_min = -1./3.; // Gev/c^(-1)
	
	
	const double t_phi = 0.005;
	const double t_deltaPhi = 0.005;
	
	const double t_eta = 1.0;
	const double t_deltaEta = 1.0;

	const double t_invPt_max = 1./3.; // Gev/c^(-1)  
	const double t_invPt_min = 0.; // Gev/c^(-1)

	
	const double t_invPt_mean = (t_invPt_max+t_invPt_min)/2.;// used for hit X correction

	const double t_x0 = 0.0;
	const double t_y0 = 0.0;
	const double t_z0 = 0.0;
	const double t_t0 = 0.0;
	const double t_deltaX0 = 0.0;
	const double t_deltaY0 = 0.0;
	const double t_deltaZ0 = 0.0;
	const double t_deltaT0 = 0.0;// mm

	 
	// Parameters for longitudinal binning (Phi0 - Pz)
	
	const double phi0Center = 0.;
	const double phi0Delta = 0.002;
	const double invPzCenter = 1/10.;
	const double invPzDelta =  0.002;
	


	// DETECTORS	

	const static unsigned nBarrels = 10;
	const static unsigned nDiscs = 30;	

	Barrel _B, B[nBarrels];
	Disc _D, D[nDiscs];

	
	
	
	Geometry() { // default constructor
	
		double xphiSigmaB = 0.0;// mm
		double xphiSigmaD = 0.0;// mm
		double zSigmaB = 0.0;// mm
		double rSigmaD = 0.0;// mm
		
		//double xphiSigmaB = 0.1;// mm
		//double xphiSigmaD = 0.1;// mm
		//double zSigmaB = 0.1;// mm
		//double rSigmaD = 0.1;// mm
		
		double tSigmaB = 0.0;
		double tSigmaD = 0.0;
		
		//double tSigmaB = 15.0;// 15 mm == 50 ps
		//double tSigmaD = 15.0;// 15 mm == 50 ps

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


		cout << " GEOMETRY INITIALIZED" << endl;
		if(iB == nBarrels && iD == nDiscs) return;
		else cout << " GEOMETRY INITIALIZATION ERROR" << endl;
		return;
		
		
	}; // end default constructor
	
	
	
	
	void print(std::ostream &out) {
	
		out << endl << "GEOMETRY DATA" << endl;
		
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


}; // end class Geometry





