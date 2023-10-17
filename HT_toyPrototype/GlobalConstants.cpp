//
//  GlobalConstants.cpp
//  MuonColliderToy
//
//  Created by Luciano Ristori on 4/28/23.
//



// GLOBAL VARIABLES

const static double Pi = 3.14159265;
const static double SpeedOfLight = 0.000299792458;// m/ps

bool Special = false;// used in debugging
	
// PARTICLE MASSES in GeV/c2

static const unsigned nMasses = 5;
double mass[] = {
	0.,
	0.000511, // electron: 1
	0.139570, // charged pion: 2
	0.493677, // charged kaon: 3
	0.938272  // proton: 4
	};
	
