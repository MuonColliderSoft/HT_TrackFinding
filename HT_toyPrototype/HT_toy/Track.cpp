//
//  Track.cpp
//  MuonColliderToy
//
//  Created by Luciano Ristori on 6/25/21
//

#include "Track.h"

extern     std::default_random_engine generator;
extern     std::uniform_real_distribution<double> distribution;
extern     std::normal_distribution<double> gauss;

//extern Parameters par;


// constructors


Track::Track(){}; // default constructor


Track::Track(TrackGeometry &g){ // constructor with parameters from specific geometry file
    
    // create random track parameters
    
    int massInd_ = 2; // it's a pion 
    //massInd_ = 3; // it's a kaon


    double x0_ = g.t_x0 + g.t_deltaX0*gauss(generator); // x at primary vertex
    double y0_ = g.t_y0 + g.t_deltaY0*gauss(generator); // y at primary vertex
    double z0_ = g.t_z0 + g.t_deltaZ0*gauss(generator); // z at primary vertex
    double t0_ = g.t_t0 + g.t_deltaT0*gauss(generator); // t at primary vertex
    
    double eta_ = 2*g.t_deltaEta * distribution(generator) + g.t_eta - g.t_deltaEta; // eta of track flat distribution
    double phi_ = 2*g.t_deltaPhi * distribution(generator) + g.t_phi - g.t_deltaPhi; // phi of track flat distribution
    
    double invPt_ =  g.t_invPt_min + (g.t_invPt_max - g.t_invPt_min) * distribution(generator); // inverse pt of track
    //if(distribution(generator) > 0.5) invPt_ = -invPt_;// choose charge// obsolete
       
    init(massInd_, x0_, y0_, z0_ , t0_, invPt_, eta_, phi_); // initialize track
    
}; // end constructor


// specific track constructor

Track::Track(int massInd_, double x0_, double y0_, double z0_, double t0_, double invPt_, double eta_, double phi_){
    
    init(massInd_, x0_, y0_, z0_ , t0_, invPt_, eta_, phi_); // initialize track
    
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
    
    
    hitList.clear();
    
    
    massInd = massInd_;
    
    // coordinates of primary vertex
    
    x0 = x0_;
    y0 = y0_;
    z0 = z0_;
    t0 = t0_;
    
    
    // primary parameters from constructor arguments
    
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
    double E2 = P2 + mass[massInd]*mass[massInd];
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


bool Track::xyDisc(DetectorGeometry &g, int iDisc, double &X, double &R, double &T){
    
    // R is the radial position of the measured hit
    // X is the distance from phi = 0 measured along the circumference
    // measurement errors are included
    
    double x,y; // cartesian coordinates of the hit
    
    // zDet is the position in z of the disc (parallel to xy)s
    
    double zDet = g.D[iDisc].z;
    
    if(cosTheta == 0.) return false; // infinite looperßß
    
    //cout << zDet << " " << z0 << " " << cosTheta << endl;
    
    double sMax = (zDet-z0)/cosTheta; // temp definition - redefined below
    
    //if(iDisc == 19) cout << zDet << " " << z0 << " " << cosTheta << endl;
    
    if (sMax <0.) return false; // track going in the wrong direction
    
    //if(iDisc == 19) cout << "24 passed" << endl;
    
    sMax = (zDet-z0)*tgTheta;
    
    if(sMax > fabs(Pi/c)) return false; // looper - ignored
    
    //T = t0 + sMax*sqrt(1.+(1./(cotTheta*cotTheta)))/beta +  g.D[iDisc].tPrec*gauss(generator);
    T = t0 + sMax*sqrt(1.+cotTheta*cotTheta)/beta +  g.D[iDisc].tPrec*gauss(generator);
    
    
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
    
    double errTang = g.D[iDisc].xphiPrec*gauss(generator); // measurement error along PHI
    double errRadial = g.D[iDisc].rPrec*gauss(generator); // measurement error along r   
    
    R = r + errRadial;
    X = R*PHI + errTang;
    

    
    // check for detector boundaries
    
    if(R > g.D[iDisc].rMax) return false;
    if(R < g.D[iDisc].rMin) return false;
       
    return true;
    
    
}; // end xyDisc

bool Track::phizBarrel(DetectorGeometry &g, int iBarrel, double &hitPhi, double &hitZ, double &hitT){
    
    // find intersection with cylindrical barrel with successive approximations using xzBarrel
    // measurement errors are included
    
    double r = g.B[iBarrel].r;
    
    if(r < sqrt(x0*x0 + y0*y0)) return false; // do not track if vertex outside barrel
    
    double rotPhi = Pi/2. - phi;
    rotate(rotPhi);
    double err = 1.;
    while(err > 1.e-8){
        //cout << " X ";
        double x_;
        if(!xzBarrel(r, x_, hitZ, hitT)){
            rotate(-rotPhi);
            return false;
        }
        else {
            double deltaPhi = atan(x_/r);
            rotate(deltaPhi);
            rotPhi += deltaPhi;
            err = fabs(deltaPhi);
            //cout << err;
        }
    }
    
    hitPhi = Pi/2. - rotPhi;
    if(hitPhi > +Pi) hitPhi -= 2*Pi;
    if(hitPhi < -Pi) hitPhi += 2*Pi;
    hitPhi *= g.B[iBarrel].r; //transform angle into linear coordinate
    
    // add measurement errors
    
    hitPhi += g.B[iBarrel].xphiPrec*gauss(generator);
    hitZ += g.B[iBarrel].zPrec*gauss(generator);
    hitT += g.B[iBarrel].tPrec*gauss(generator); 
    
    rotate(-rotPhi);
    
    // check for z boundaries
    
    if(hitZ > g.B[iBarrel].zMax) return false;
    if(hitZ < g.B[iBarrel].zMin) return false;
   
       
    return true;
};


