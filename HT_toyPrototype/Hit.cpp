//
//  Hit.cpp
//  MuonColliderToy
//
//  Created by Luciano Ristori on 2/28/23
//


        
        Hit::Hit(char _hitType, int _iLayer, double _x1, double _x2, double _t, unsigned _trackInd){
		
			hitType = _hitType;
			iLayer = _iLayer;
			x1 = _x1;
			x2 = _x2;
			t = _t;
			trackInd = _trackInd;
			layerInd = iLayer;
			if(hitType == 'D')layerInd = iLayer + dg.nBarrels;
        
        }
            
		void Hit::print(ostream &out){
	
		out <<  "type " << hitType << " layer: " << layerInd  << " x1: " << x1 << " x2: " << x2 << " time: " << t << " iT: " <<  trackInd << std::endl;
	
		}
		
		void Hit::write(ostream &out){
		
		double X, Y, Z;
			
		XYZ(dg, X,Y,Z);
	
		out << hitType << " " << layerInd  << " " << X << " " << Y << " " << Z << " " << t << " " <<  trackInd << std::endl;
	
		}
		
		void Hit::XYZ(DetectorGeometry &g, double &X, double &Y, double &Z){
			double PHI;
			if(hitType == 'B'){
				Z = x2;
				double R = g.B[iLayer].r;
				double PHI = x1/R;
				X = R*cos(PHI);
				Y = R*sin(PHI);	
			}
			else {
				Z = g.D[iLayer].z;
				PHI = x1/x2;
				X = x2*cos(PHI);
				Y = x2*sin(PHI);		
			}
			
		}
		
		double Hit::realPhi(DetectorGeometry &g){
			double PHI;
			if(hitType == 'B') {
				double R = g.B[iLayer].r;
				PHI = x1/R;
			}
			else PHI = x1/x2;
			
			return PHI;			
		}
		
    
		void Hit::printXYZ(DetectorGeometry &g, ostream &out){
		
			double X, Y, Z;
			
			XYZ(g, X,Y,Z);
	
			out << "{" << X << "," << Y << "," << Z << "}";
		}
		
		
		double Hit::timeExpected(DetectorGeometry &g, double mass, double invPt){
		
		// TOF of a particle of mass "mass" and inverse pt "invPt" from the origin to the hit point
		// This takes into account both the curvature of the track and the beta.
		// r and z have the usual meaning in the global coordinate system
		
			const static double k= 3333.333333333/par.magneticField;; // mm/GeV  (Radius/Pt)
			
			double r = 0.;
			double z = 0.;
			
			
			if(hitType == 'B'){
				z = x2;
				r = g.B[iLayer].r;	
			}
			
			if(hitType == 'D'){
				z = g.D[iLayer].z;
				r = x2;	
			}
			
			return fabs(2*k*sqrt(mass*mass + 1./invPt/invPt*(1 + z*z/r/r))*asin(invPt*r/(2.*k)));
		
		}
    
    	bool Hit::isSeed(){
		
			if(layerInd >= 4 && layerInd <= 9) return true;
			if(layerInd >= 14 && layerInd <= 24) return true;
			if(layerInd >= 29 && layerInd <= 39) return true;
			return false;
        
        }
    
    