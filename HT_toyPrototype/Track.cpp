//
//  Track.cpp
//  MuonColliderToy
//
//  Created by Luciano Ristori on 6/25/21
//

#include "Track.h"

extern double phi_to_xi(double phi);
extern double eta_to_xj(double eta);
extern double invPt_to_xk(double invPt);

// constructors


Track::Track(){}; // default constructor


Track::Track(TrackGeometry &g){ // constructor with parameters from specific geometry file
    
    // create random track parameters
    
     int massInd_ = 2;// it's a pion 
     
   // if(distribution(generator_trk) > 0.5) massInd_ = 2; // it's a pion 
    	//else massInd_ = 3; // it's a kaon
     
		//else massInd_ = 4; // it's a proton

    double x0_ = g.t_x0 + g.t_deltaX0*gauss(generator_trk); // x at primary vertex
    double y0_ = g.t_y0 + g.t_deltaY0*gauss(generator_trk); // y at primary vertex
    double z0_ = g.t_z0 + g.t_deltaZ0*gauss(generator_trk); // z at primary vertex
    double t0_ = g.t_t0 + g.t_deltaT0*gauss(generator_trk); // t at primary vertex
    
    double eta_ = 2*g.t_deltaEta * distribution(generator_trk) + g.t_eta - g.t_deltaEta; // eta of track flat distribution
    double phi_ = 2*g.t_deltaPhi * distribution(generator_trk) + g.t_phi - g.t_deltaPhi; // phi of track flat distribution
    
    double invPt_ =  g.t_invPt_min + (g.t_invPt_max - g.t_invPt_min) * distribution(generator_trk); // inverse pt of track
       
    init(massInd_, x0_, y0_, z0_ , t0_, invPt_, eta_, phi_); // initialize track
    
}; // end constructor


// specific track constructors

Track::Track(int massInd_, double x0_, double y0_, double z0_, double t0_, double invPt_, double eta_, double phi_){
    
    init(massInd_, x0_, y0_, z0_ , t0_, invPt_, eta_, phi_); // initialize track
    
};

Track::Track(double mass_, double x0_, double y0_, double z0_, double t0_, double invPt_, double eta_, double phi_){
    
    init(mass_, x0_, y0_, z0_ , t0_, invPt_, eta_, phi_); // initialize track
    
};



void Track::print(std::ostream &out, int mode = 0){
	// mode = 0: track parameters only
	// mode = 1: list all hits
    
    out << "ID:" << ID << " ind: " << massInd << " beta: " << beta << " vertex: x: "  
        << setw(8)<< x0  << " y: "  << setw(8)<< y0  << " z: "  << setw(8)<< z0 << " t: " << setw(8) << t0
        << " charge: " << charge << " pt: "  << setw(12)<< Pt  << " eta: " << setw(12) << eta 
        << " phi: "  << setw(12)<< phi << endl;
        
    if(mode == 0) return;
    
    for(unsigned iH = 0; iH != hitList.size(); ++iH)  (hitList[iH]).print(out);
    
}; // end print



// common code

void Track::init(int massInd_, double x0_, double y0_, double z0_, double t0_, double invPt_, double eta_, double phi_){
	init(Mass[massInd_], x0_, y0_, z0_, t0_, invPt_, eta_, phi_);
};


	
void Track::init(double mass_, double x0_, double y0_, double z0_, double t0_, double invPt_, double eta_, double phi_){
      
    hitList.clear();   
    
   // massInd = massInd_;
    
    // coordinates of primary vertex
    
    x0 = x0_;
    y0 = y0_;
    z0 = z0_;
    t0 = t0_;
    
    
    // primary parameters from constructor arguments
    
    mass = mass_;
    invPt = invPt_;
    eta = eta_;
    phi = phi_; 
    
    // derived parameters
    
    c = 3.e-4*invPt*par.magneticField; // signed curvature in mmm^(-1)
    
    charge = 1.;
    	if(invPt < 0.) charge = -1.;
    
    cotTheta = sinh(eta); // cotTheta = pz/pt
    tgTheta = 1./cotTheta;
    
    Pt = 1./invPt; // transverse momentum
    	if(Pt < 0.) Pt = -Pt;
    Pz = Pt*cotTheta; // longitudinal momentum
    double P2 = Pt*Pt + Pz*Pz;
    double E2 = P2 + mass*mass;
    beta = sqrt(P2/E2);// velocity
    
    cosTheta = Pz/sqrt(P2);
    
    phi0 = phi - c*Pt*z0/2./Pz;
    
    
}; // end init


// rigid rotation of track around z axis

void Track::rotate(double phiRot){
    
    double xTemp = x0*cos(phiRot) - y0*sin(phiRot);
    double yTemp = x0*sin(phiRot) + y0*cos(phiRot);
    x0 = xTemp;
    y0 = yTemp;
    phi += phiRot;
    if(phi > +Pi) phi -= 2*Pi;
    if(phi < -Pi) phi += 2*Pi;
    
    return;
    
}; // end rotate


// calculate exact coordinates of hits for flat barrel and disc layers
// return "false" if no intersection
// measurement errors are NOT included

bool Track::xzBarrel(double yDet, double &x, double &z, double &t){
    
    // yDet is the position of the barrel plane (parallel to xz)
    
    double arg = c*(y0-yDet) + cos(phi);
    
    //cout << "*** c: " << c << " y0: " << y0 << " yDet: " << yDet << endl;
    
    if(fabs(arg) > 1.)return false;
    
    if(fabs(c) > cMin){
        
        // find two solutions for s (range is (-2*Pi/c, +2*Pi/c)
        // sol is the arc length measured on the transverse plane
        // starting from the origin of the track at y0
        
        double sol1 = (+acos(arg) - phi)/c;
        double sol2 = (-acos(arg) - phi)/c;
        
        // bring solutions to positive s
        
        if(sol1 < 0.) sol1 += fabs(2*Pi/c);
        if(sol2 < 0.) sol2 += fabs(2*Pi/c);
        
        double sMax = min(sol1,sol2);
        //cout << "*** sol1: " << sol1 << " sol2: " << sol2 << endl;
        
        x = x0 + 1./c*(sin(c*sMax + phi)-sin(phi));
        z = z0 + cotTheta*sMax;
        t = t0 + sMax*sqrt(1.+cotTheta*cotTheta)/beta;
        return true;
    }
    else {
        // very large momentum
        // use first order approximation in c (curvature)
        
        if(fabs(sin(phi)) < 0.001) return false; // no intersection with barrel
        
        double sMax = (yDet - y0)/sin(phi) - (yDet - y0)*(yDet - y0)*cos(phi)/2./sin(phi)/sin(phi)/sin(phi)*c;
        x = x0 + sMax*cos(phi) - 0.5*sMax*sMax*c*sin(phi);
        z = z0 + cotTheta*sMax;
        t = t0 + sMax*sqrt(1.+cotTheta*cotTheta)/beta;
        
        return true;
    }
    
};


bool Track::xyDisc(DetectorGeometry &g, int iDisc, double &X, double &R, double &T, bool smear = false, bool checkBoundaries = true){
    
    // R is the radial position of the measured hit
    // X is the distance from phi = 0 measured along the circumference
    // measurement errors are included
    
    double x,y; // cartesian coordinates of the hit
    
    // zDet is the position in z of the disc (parallel to xy)s
    
    double zDet = g.D[iDisc].z;
    
    if(cosTheta == 0.) {
    	cout << "infinite looper" << endl;
    	return false; // infinite looper   	
    }
    
    //cout << zDet << " " << z0 << " " << cosTheta << endl;
    
    double sMax = (zDet-z0)/cosTheta; // temp definition - redefined below
    
    if (sMax <0.) {
    	//cout << "track going in wrong direction iDisc " << iDisc << endl;
    	return false; // track going in the wrong direction
    }
    
    //if(iDisc == 19) cout << "24 passed" << endl;
    
    sMax = (zDet-z0)*tgTheta;
    
    if(sMax > 2*fabs(Pi/c)) {
    	//cout << "looper" << endl;
    	return false; // looper - ignored   	
    }
    	
    
    //T = t0 + sMax*sqrt(1.+(1./(cotTheta*cotTheta)))/beta +  g.D[iDisc].tPrec*gauss(generator_trk);
    T = t0 + sMax*sqrt(1.+cotTheta*cotTheta)/beta; 
    
    
    if(fabs(c) > cMin){
        x = x0 + 1./c*(sin(c*sMax + phi) - sin(phi));
        y = y0 - 1./c*(cos(c*sMax + phi) - cos(phi));
    }
    else {
        // very large momentum
        // use first order approximation in c (curvature)
        
        x = x0 + sMax*cos(phi) - 0.5*sMax*sMax*c*sin(phi);
        y = y0 + sMax*sin(phi) + 0.5*sMax*sMax*c*cos(phi);
    }
    
    
    // find polar coordinates
    
    double r = sqrt(x*x + y*y);
    double PHI = atan2(y,x);
    
    if(smear){
     
		double errTang = g.D[iDisc].xphiPrec*gauss(generator_trk); // measurement error along PHI
		double errRadial = g.D[iDisc].rPrec*gauss(generator_trk); // measurement error along r   
	
		R = r + errRadial;
		X = R*PHI + errTang;
		T += g.D[iDisc].tPrec*gauss(generator_trk);
    }
    else {
		R = r;
		X = R*PHI;
    }

	if(!checkBoundaries) return true; // 
    
    // check for detector boundaries
    
    if(R > g.D[iDisc].rMax) {
    	//cout << "Disc " << iDisc << "out of high R bounds" << endl;
    	return false;
    }
    if(R < g.D[iDisc].rMin) {
    	//cout << "Disc " << iDisc << "out of low R bounds" << endl;
    	return false;
    }
       
    return true;
    
    
}; // end xyDisc

bool Track::phizBarrel(DetectorGeometry &g, int iBarrel, double &hitPhi, double &hitZ, double &hitT, bool smear = false, bool checkBoundaries = true){
    
    // find intersection with cylindrical barrel with successive approximations using xzBarrel
    // measurement errors are included
    
    double r = g.B[iBarrel].r;
    
    if(r < sqrt(x0*x0 + y0*y0)) return false; // do not track if vertex outside barrel
    
    double rotPhi = Pi/2. - phi;
    rotate(rotPhi);
    double err = 1.;
    while(err > 1.e-8){
        double x_;
        if(!xzBarrel(r, x_, hitZ, hitT)){
            rotate(-rotPhi);
            cout << "xzBarrel fails. rotPhi = " << rotPhi << endl;
            return false;
        }
        else {
            double deltaPhi = atan(x_/r);
            rotate(deltaPhi);
            rotPhi += deltaPhi;
            err = fabs(deltaPhi);
        }
    }
    
    hitPhi = Pi/2. - rotPhi;
    if(hitPhi > +Pi) hitPhi -= 2*Pi;
    if(hitPhi < -Pi) hitPhi += 2*Pi;
    hitPhi *= g.B[iBarrel].r; //transform angle into linear coordinate
    
    // add measurement errors
        
    if(smear){
		hitPhi += g.B[iBarrel].xphiPrec*gauss(generator_trk);
		hitZ += g.B[iBarrel].zPrec*gauss(generator_trk);
		hitT += g.B[iBarrel].tPrec*gauss(generator_trk); 
    }
    
    rotate(-rotPhi);
    
    if(!checkBoundaries) return true;
    
    // check for z boundaries
    
    if(hitZ > g.B[iBarrel].zMax) {
    	//cout << " Barrel "  << iBarrel << " Zmax limit exceeded" << endl;
    	return false;
    	}
    if(hitZ < g.B[iBarrel].zMin) {
    	//cout << " Barrel "  << iBarrel << " Zmin limit exceeded" << endl;
    	return false;
    	}
   
       
    return true;
    
}; 


   

double Track::hitChi2(DetectorGeometry &g, Hit &h){

	double x1, x2, t;
	double x1err, x2err, terr;
	unsigned iLayer = h.iLayer;
	bool noSmear = false;
	bool checkBoundaries = false;
	
	if(h.hitType == 'B'){
		if(!phizBarrel(g, h.iLayer, x1, x2, t, noSmear, checkBoundaries)){
			cout << "phizBarrel fails. Barrel " << iLayer << endl;
			h.print(cout);
		}
		x1err = g.B[iLayer].xphiPrec;
		x2err = g.B[iLayer].zPrec;
		terr = g.B[iLayer].tPrec;		
	}
	else {
		if(!xyDisc(g, h.iLayer, x1, x2, t, noSmear, checkBoundaries)) {
			cout << "xyDisc fails. Disc " << iLayer << endl; 
			h.print(cout);
		}
		x1err = g.D[iLayer].xphiPrec;
		x2err = g.D[iLayer].rPrec;
		terr = g.D[iLayer].tPrec;
	}
		
	//cout << "*** Track::hitChi2 ***" << endl;
	//h.print(cout);
/*	print(cout);
	cout << " iLayer: " << iLayer << endl;
	cout                  <<  "track: " << x1   << "	" << x2   <<"	"<< t   << endl;
	cout << " " << h.hitType << "hit: " << h.x1 << "	" << h.x2 <<"	"<< h.t << endl;
	*///cout << "**********************"<< endl;
	
	double  chi2 =  (x1-h.x1)/x1err*(x1-h.x1)/x1err; 
			chi2 += (x2-h.x2)/x2err*(x2-h.x2)/x2err; 
			chi2 += (t-h.t)/terr*(t-h.t)/terr;
    return chi2; 
};

void Track::getIJK(unsigned &I, unsigned &J, unsigned &K){

	I = phi_to_xi(phi);
	J = eta_to_xj(eta);
	K = invPt_to_xk(invPt);
	
	return;

}; // returns indices in HTM array


