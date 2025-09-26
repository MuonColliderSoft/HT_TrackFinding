
#ifndef MDHT_HIT_CPP
#define MDHT_HIT_CPP

//
//  Hit.cpp
//  MuonColliderToy
//
//  Created by Luciano Ristori on 2/28/23
//  Modified for use with Muon Collider code April 3 2025
//

#include "Hit.h"

//////////////////////////////////////////////////////////////////////////////////////     
            
        // dummy hit constructor
       
        Hit::Hit(char _hitType){
        	hitType = _hitType;
        }
        
////////////////////////////////////////////////////////////////////////////////////// 

		// Regular Hit constructor  
		  
        Hit::Hit(char _hitType, int _CellID, int _layerInd, int _trackInd, 
        	double _x1, double _x2, double _t, double _x, double _y, double _z,
        	double _u1x, double _u1y, double _u1z, double _u2x, double _u2y, double _u2z,
        	double _errx1, double _errx2){
        
			hitType = _hitType;
			CellID = _CellID;
			layerInd = _layerInd;
			trackInd = _trackInd;
			
			x1 = _x1;
			x2 = _x2;
			t = _t; 
			x = _x;
			y = _y;
			z = _z; 
			
			u1x = _u1x;
			u1y = _u1y;
			u1z = _u1z;
			
			u2x = _u2x;
			u2y = _u2y;
			u2z = _u2z;
			
			errx1 = _errx1;
			errx2 = _errx2;
						   
        }
        
////////////////////////////////////////////////////////////////////////////////////// 
  
		void Hit::print(ostream &out){
	
		out << " ID " << ID << " type: " << hitType << " CellID: " << CellID 
			<< " layerInd: " << layerInd  << " trackInd: " << trackInd 
			<< " x1: " << x1 << " x2: " << x2 << " time: " << t  
			<< " x: " << x  << " y: " << y  << " z: " << z  
			<< " u1x: " <<  u1x << " u1y: " <<  u1y << " u1z: " <<  u1z 
			<< " u2x: " <<  u2x << " u2y: " <<  u2y << " u2z: " <<  u2z 
			<< " errx1: " << errx1 << " errx2: " << errx2 << std::endl;		
	
		}
		
////////////////////////////////////////////////////////////////////////////////////// 
		
		void Hit::write(ostream &out){
	
			out << hitType << " " << CellID << " " << layerInd  << " " << trackInd << " " 
				<< x1 << " " << x2 << " " << t  << " " 
				<< x  << " " << y  << " " << z  << std::endl;
	
		}
			
////////////////////////////////////////////////////////////////////////////////////// 
    
    	bool Hit::isSeed(){
    		return true;

    		const unsigned int system = (unsigned) ( CellID & 0x1f ); 		 
    		if(system >= 3) return true;
			return false;
        
        }
        
////////////////////////////////////////////////////////////////////////////////////// 


#endif //  MDHT_HIT_CPP    
    