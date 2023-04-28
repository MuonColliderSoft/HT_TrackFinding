//
//  TrackGeometry.cpp
//  MuonColliderToy
//
//  Created by Luciano Ristori on 4/28/23.
//


class TrackGeometry {


	public:
			
    
    // TRACK GENERATION PARAMETERS
  
	
	double t_phi;
	double t_deltaPhi;
	
	double t_eta;
	double t_deltaEta;

	double t_invPt_max;  
	double t_invPt_min;

	
	double t_x0;
	double t_y0;
	double t_z0;
	double t_t0;
	double t_deltaX0;
	double t_deltaY0;
	double t_deltaZ0;
	double t_deltaT0;
	
	TrackGeometry() { // default constructor
	
		t_phi = par.geo_def_t_phi;
		t_deltaPhi = par.geo_def_t_deltaPhi;
	
		t_eta = par.geo_def_t_eta;
		t_deltaEta = par.geo_def_t_deltaEta;

		t_invPt_max = par.geo_def_t_invPt_max;  
		t_invPt_min = par.geo_def_t_invPt_min;

	
		t_x0 = par.geo_def_t_x0;
		t_y0 = par.geo_def_t_y0;
		t_z0 = par.geo_def_t_z0;
		t_t0 = par.geo_def_t_t0;
		t_deltaX0 = par.geo_def_t_deltaX0;
		t_deltaY0 = par.geo_def_t_deltaY0;
		t_deltaZ0 = par.geo_def_t_deltaZ0;
		t_deltaT0 = par.geo_def_t_deltaT0;	
	
	}



}; // end class TrackGeometry





