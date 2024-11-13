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
		_B.r = 1510.; _B.zMax = 1360.; _B.zMin =-_B.zMax; _B.xphiPrec = xphiSigmaB; _B.zPrec = zSigmaB; _B.tPrec = tSigmaB;    B[iB++] = _B;
		
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


//Resolutions from Multiple Scattering

if(true){
	
	B[0].xphiPrec = 0.01;
	B[0].zPrec = 0.1;
	B[1].xphiPrec = 0.0148578;
	B[1].zPrec = 0.100793;
	B[2].xphiPrec = 0.0271705;
	B[2].zPrec = 0.103492;
	B[3].xphiPrec = 0.0474344;
	B[3].zPrec = 0.110331;
	B[4].xphiPrec = 0.0520281;
	B[4].zPrec = 0.120851;
	B[5].xphiPrec = 0.306575;
	B[5].zPrec = 0.356449;
	B[6].xphiPrec = 0.559915;
	B[6].zPrec = 0.605871;
	B[7].xphiPrec = 0.916051;
	B[7].zPrec = 0.990146;
	B[8].xphiPrec = 1.54581;
	B[8].zPrec = 1.44297;
	B[9].xphiPrec = 2.34229;
	B[9].zPrec = 1.95837;
	B[10].xphiPrec = 0.01;
	B[10].zPrec = 0.1;
	D[0].xphiPrec = 0.392332;
	D[0].rPrec = 0.563409;
	D[1].xphiPrec = 0.17323;
	D[1].rPrec = 0.234711;
	D[2].xphiPrec = 0.0398203;
	D[2].rPrec = 0.108302;
	D[3].xphiPrec = 0.030051;
	D[3].rPrec = 0.104338;
	D[4].xphiPrec = 0.1826;
	D[4].rPrec = 0.233638;
	D[5].xphiPrec = 0.25901;
	D[5].rPrec = 0.304069;
	D[6].xphiPrec = 0.291457;
	D[6].rPrec = 0.32612;
	D[7].xphiPrec = 0.293327;
	D[7].rPrec = 0.322848;
	D[8].xphiPrec = 0.300988;
	D[8].rPrec = 0.327248;
	D[9].xphiPrec = 0.299469;
	D[9].rPrec = 0.323565;
	D[10].xphiPrec = 0.289007;
	D[10].rPrec = 0.311966;
	D[11].xphiPrec = 1.12948;
	D[11].rPrec = 1.41284;
	D[12].xphiPrec = 1.14543;
	D[12].rPrec = 1.34727;
	D[13].xphiPrec = 1.14881;
	D[13].rPrec = 1.30031;
	D[14].xphiPrec = 1.17838;
	D[14].rPrec = 1.29757;
	D[15].xphiPrec = 0.390393;
	D[15].rPrec = 0.561005;
	D[16].xphiPrec = 0.171872;
	D[16].rPrec = 0.233055;
	D[17].xphiPrec = 0.0398903;
	D[17].rPrec = 0.10832;
	D[18].xphiPrec = 0.0304714;
	D[18].rPrec = 0.104477;
	D[19].xphiPrec = 0.182193;
	D[19].rPrec = 0.233169;
	D[20].xphiPrec = 0.259053;
	D[20].rPrec = 0.304113;
	D[21].xphiPrec = 0.29166;
	D[21].rPrec = 0.326195;
	D[22].xphiPrec = 0.29714;
	D[22].rPrec = 0.326774;
	D[23].xphiPrec = 0.302233;
	D[23].rPrec = 0.328492;
	D[24].xphiPrec = 0.299268;
	D[24].rPrec = 0.32335;
	D[25].xphiPrec = 0.290106;
	D[25].rPrec = 0.313087;
	D[26].xphiPrec = 1.12677;
	D[26].rPrec = 1.40861;
	D[27].xphiPrec = 1.14634;
	D[27].rPrec = 1.34849;
	D[28].xphiPrec = 1.14631;
	D[28].rPrec = 1.2972;
	D[29].xphiPrec = 1.17679;
	D[29].rPrec = 1.29571;
}

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





