//
//  HTArray.cpp
//  MuonColliderToy
//
//  Created by Luciano Ristori on 5/23/22
//


 extern int TrainingPhase;
 extern bool Special;

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

    class Hitstat {
    
    	// This is the class that contains all that is needed for each layer in each HT cell
    
    	    
		public:
		
		bool hitLayer; // this layer has been hit in this HT element	
		vector<unsigned> hitIDList; // list of hit ID's (indices in hit vector in event)
    
    	long int nEntries;
    	double mean[3]; // mean values of x1, x2 and t
    	double low[3]; // min value
    	double high[3]; // max value
    	double c00, c11, c22, c01, c02, c12; //  covariance matrix c[i,j]
    	double rotAngle; // rotation angle for x2-t diagonalization
    	double sinAngle; // sin of rotAngle
    	double cosAngle; // cos of rotAngle
    	double u0, u1, u2; // diagonal coordinates
    	double u0m, u1m, u2m; // mean of diagonal coordinates
    	double u0v, u1v, u2v; // variance of diagonal coordinates
    	double u0Max, u0Min, u1Max, u1Min, u2Max, u2Min; // min-max values
    	
    	
    	Hitstat(){ // constructor
    		nEntries = 0;
    		hitIDList.clear();
    		c00 = 0.;
    		c11 = 0.;
    		c22 = 0.;
    		c01 = 0.;
    		c02 = 0.;
    		c12 = 0.;
    		for(int k = 0; k!= 3; ++k) {
    			mean[k] = 0.;
    			low[k] = 0.;
    			high[k] = 0.;
    		}
    		rotAngle = 0.;
    		sinAngle = sin(rotAngle);
    		cosAngle = cos(rotAngle);
    		u0 = 0.;
    		u1 = 0.;
    		u2 = 0.;
    		u0m = 0.;
    		u1m = 0.;
    		u2m = 0.;
    		u0v = 0.;
    		u1v = 0.;
    		u2v = 0.;
    		u0Max = 0.;
    		u0Min = 0.;
    		u1Max = 0.;
    		u1Min = 0.;
    		u2Max = 0.;
    		u2Min = 0.;
    	}
    
    	void train(Hit &h){
    	
    	
    	// incremental calculation of means, variances and covariance.
    	// classic algorithm
    	
   			++ nEntries;
   			
			if(TrainingPhase == 1){
			
				mean[0] += (h.x1 - mean[0])/nEntries;
				mean[1] += (h.x2 - mean[1])/nEntries;
				mean[2] += (h.t  - mean[2])/nEntries;
			
				if(nEntries > 1 ){
					c00 += (h.x1 - mean[0])*(h.x1 - mean[0])/(nEntries-1) - c00/nEntries;
					c11 += (h.x2 - mean[1])*(h.x2 - mean[1])/(nEntries-1) - c11/nEntries;
					c22 += (h.t  - mean[2])*(h.t  - mean[2])/(nEntries-1) - c22/nEntries;
					c01 += (h.x1 - mean[0])*(h.x2 - mean[1])/(nEntries-1) - c01/nEntries;
					c02 += (h.x1 - mean[0])*(h.t  - mean[2])/(nEntries-1) - c02/nEntries;
					c12 += (h.x2 - mean[1])*(h.t  - mean[2])/(nEntries-1) - c12/nEntries;
				}
				
				if(nEntries == 1 ){
					low[0] = h.x1;
					high[0] = h.x1;
					low[1] = h.x2;
					high[1] = h.x2;
					low[2] = h.t;
					high[2] = h.t;
				} else {
					if(h.x1 > high[0]) high[0] = h.x1;
					if(h.x1 < low[0]) low[0] = h.x1;
					if(h.x2 > high[1]) high[1] = h.x2;
					if(h.x2 < low[1]) low[1] = h.x2;
					if(h.t > high[2]) high[2] = h.t;
					if(h.t < low[2]) low[2] = h.t;	
				}
				
				
			} // end if(TrainingPhase == 1)
			
			
			
			if(TrainingPhase >= 2){
				// rotation into diagonal space u0,u1,u2
			
				u0 = h.x1 - mean[0];
				u1 = (h.x2-mean[1]) * cosAngle - (h.t-mean[2]) * sinAngle;
				u2 = (h.x2-mean[1]) * sinAngle + (h.t-mean[2]) * cosAngle;
			
				// copy vector u into hit
			
				h.u0 = u0;
				h.u1 = u1;
				h.u2 = u2;
				
			}// end if(TrainingPhase >= 2)
				
			if(TrainingPhase == 2) {
			
				// find min and max values of u's
			 
				if(nEntries == 1){
					u0Max = u0;
					u0Min = u0;
					u1Max = u1;
					u1Min = u1;
					u2Max = u2;
					u2Min = u2;
				} else {
					if(u0 > u0Max) u0Max = u0;
					if(u0 < u0Min) u0Min = u0;
					if(u1 > u1Max) u1Max = u1;
					if(u1 < u1Min) u1Min = u1;
					if(u2 > u2Max) u2Max = u2;
					if(u2 < u2Min) u2Min = u2;	
				}
			
				// incremental calculation of mean and variance of diagonal coordinates u0,u1,u2
			
				u0m += (u0 - u0m)/nEntries;
				u1m += (u1 - u1m)/nEntries;
				u2m += (u2 - u2m)/nEntries;
			
				if(nEntries > 1 ){
					u0v += (u0 - u0m)*(u0 - u0m)/(nEntries-1) - u0v/nEntries;
					u1v += (u1 - u1m)*(u1 - u1m)/(nEntries-1) - u1v/nEntries;
					u2v += (u2 - u2m)*(u2 - u2m)/(nEntries-1) - u2v/nEntries;			
				}
			} // end if(TrainingPhase == 2)
			
   		} // end train
   		
   		void diagonalize(){
   		
   			// finds diagonalizing rotation angle
   			// for this cell, this layer
   		
   			rotAngle = -0.5*atan2(2*c12,c11-c22);
   			
   		}
    	
    	    
    	long int GetEntries(){
    		
    		return nEntries;
    	
    	}
    	
    	double GetMean(int iPar){
    	
    		return mean[iPar];
    	
    	}
    	
    	void reset(){
    	
    		hitLayer = false;
    		hitIDList.clear();
    	}
    	
    	void print(ostream &out){
    		
    		out << "N: " << nEntries << " x1: " << mean[0] << " x2: " <<  mean[1] << " t: " <<  mean[2] << " ";
    		//out << "c00: " << c00 << " c11: " << c11 << " c22: " << c22 << " c01: " << c01 << " c02: " << c02 << " c12: " << c12; 
    		out << " rotAngle: " << rotAngle; 
    		out << "     u0:[" << u0Min << "," << u0Max << "]; u1:[" << u1Min << "," << u1Max << "]; u2:[" << u2Min << "," << u2Max << "]";
    	}
    	
    	void printHits(ostream &out){
    	
    		int nHits = hitIDList.size();
    		out << "nHits = " << nHits << ": ";
    		for(int iH = 0; iH != nHits; ++ iH) out << hitIDList[iH] << " ";
    		out << endl;
    	}
    	
    	
    	void writeHits(ostream &out, vector<Hit> & hitList){
    	
    		int nHits = hitIDList.size();
    		if(nHits) out << nHits << endl;
    		for(int iH = 0; iH != nHits; ++ iH) {
    			hitList[hitIDList[iH]].write(out);
    		}
    		//out << endl;
    	}
    	
    	void write(ostream &out){
    	
    		out << mean[0] << " " << mean[1] << " " << mean[2] << " ";
    		out << low[0] << " " << low[1] << " " << low[2] << " ";
    		out << high[0] << " " << high[1] << " " << high[2] << " ";
    		out << rotAngle << " " << u0Min << " " << u1Min << " " << u2Min << " " << u0Max << " " << u1Max << " " << u2Max << " ";
    	}
    	
    }; // end Hitstat
   
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

    class HTArrayElement {   
          
    // This is one element of the HT Array with parameters and accessories
           
    public:
               
        unsigned nBarrels = dg.nBarrels;       
        unsigned nHitLayers	; // number of layers hit for current event in this element
             
        // list of layers for this element implemented as a map
        map<int, Hitstat> layerIndHitStat;
        map<int, Hitstat>::iterator it;
        
        // next two params calculated in training phase and loaded at initialization
        unsigned minLayers = 10000; // min number of layers hit by a track from this param space
        unsigned maxLayers = 0;// max number of layers hit by a track from this param space
        
       // train this array element
        // add this layer to the list for this element if missing
        
        void train(Hit &h){
        	
        	layerIndHitStat[h.layerInd].train(h);
        }
        
        // diagonalize this element
        
        void diagonalize (){
        
        	for(it = layerIndHitStat.begin(); it != layerIndHitStat.end(); ++it){      
        			(it->second).diagonalize(); 
        	}      
        }
        
        // reset all hit flags in all layers of this element
        
         void reset (){
         
         	nHitLayers = 0;
        
        	for(it = layerIndHitStat.begin(); it != layerIndHitStat.end(); ++it){      
        			(it->second).reset(); 
        	}      
        }
        
        
        
        void print(ostream &out){
        	out << "min,max Layers: " << minLayers << "," << maxLayers << endl;
        	for(it = layerIndHitStat.begin(); it != layerIndHitStat.end(); ++it){
        		out << "   lay: " << it->first; out << " ";
        			(it->second).print(out); out<< endl;
        	}
        }
        
        
        void printHits(ostream &out){
        	out << "min,max Layers: " << minLayers << "," << maxLayers << endl;
        	for(it = layerIndHitStat.begin(); it != layerIndHitStat.end(); ++it){
        		out << "lay: " << it->first; out << " ";
        			(it->second).printHits(out); out << endl;
        	}
        }
        
          void writeHits(ostream &out, vector<Hit> & hitList){
        	//out << "min,max Layers: " << minLayers << "," << maxLayers << endl;
        	for(it = layerIndHitStat.begin(); it != layerIndHitStat.end(); ++it){
        		//out << "lay: " << it->first; out << " ";
        			(it->second).writeHits(out, hitList);// out << endl;
        	}
        }
        
        
             
        int fill(Hit &h, int mode){
        
        
        	
		
			it = layerIndHitStat.find(h.layerInd);
			if(it != layerIndHitStat.end()){
				Hitstat stat = it->second;
									
				double u0 =  h.x1 - stat.mean[0];
				if(Special) cout << u0 << " " << stat.u0Min << " " 	<< stat.u0Max << endl;	
				if(u0 < stat.u0Min || u0 > stat.u0Max) return -1;
					if(Special)cout << "good u0" << endl;
				double u1 = (h.x2 - stat.mean[1]) * stat.cosAngle - (h.t - stat.mean[2]) * stat.sinAngle;
				if(Special) cout << u1 << " " << stat.u1Min << " " 	<< stat.u2Max << endl;		
				if(u1 < stat.u1Min || u1 > stat.u1Max) return -2;	
					if(Special)cout << "good u1"<< endl;
				double u2 = (h.x2 - stat.mean[1]) * stat.sinAngle + (h.t - stat.mean[2]) * stat.cosAngle;	
				if(Special) cout << u2 << " " << stat.u2Min << " " 	<< stat.u2Max << endl;		
				if(u2 < stat.u2Min || u2 > stat.u2Max) return -3;	
					if(Special)cout << "good u2"<< endl;
						
				if(!stat.hitLayer)++nHitLayers;
				
				return 0;
					
			} 
			
			else return -4; 
			  
        }// end fill
        
        
        
    }; // end HTArrayElement
    
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
    
    class HTArray {
               
    // This is the whole Hough Transform Array with parameters and accessories
  
    	public:
    	
    	ofstream outfile;	
    	
    	// one histogram of a single HT cell content as an exampl
			
		TH3I *H3D_HTcellx = new TH3I("H3D_HTcellx","H3D_HTcellx",40,0.,1.,40,0.,1.,40,0.,1.);
		TH3I *H3D_HTcellu = new TH3I("H3D_HTcellu","H3D_HTcellu",40,0.,1.,40,0.,1.,40,0.,1.);

    
    	// dimensions of the array
    	
                const static unsigned NphiBins = HTA_NphiBins;
		const static unsigned NetaBins = HTA_NetaBins;
		const static unsigned NinvptBins = HTA_NinvptBins;
		
		// origin and steps of the array
		
		double phiMin;
		double phiStep;
		double etaMin;
		double etaStep;
		double invPtMin;
		double invPtStep;
		
		// HT Array bin and layer to plot content in 3D
		
		int plotBinX = par.HTA_plotBinX;
		int plotBinY = par.HTA_plotBinY;
		int plotBinZ = par.HTA_plotBinZ;
		int plotLay = par.HTA_plotLay;
		
		// HT Array proper	

                HTArrayElement ArrElem[HTA_NphiBins][HTA_NetaBins][HTA_NinvptBins];
      
		HTArray(){ // constructor
		
			phiMin = tg.t_phi - tg.t_deltaPhi;
			phiStep = 2*tg.t_deltaPhi/double(NphiBins);
			etaMin = tg.t_eta - tg.t_deltaEta;
			etaStep = 2*tg.t_deltaEta/double(NetaBins);
			invPtMin = tg.t_invPt_min;
			invPtStep = (tg.t_invPt_max - tg.t_invPt_min)/double(NinvptBins);
			
			string mapDataFileName = par.HTA_mapDataFileName;
			outfile.open(mapDataFileName);
				
		}
		
				
		void print(ostream &out, int level = 0){
		
			long int totalNbins = NphiBins*NetaBins*NinvptBins;
			out << "Total number of bins in the array: " << totalNbins << endl;
		
			out << "phi bins: " << NphiBins << "; phi min: " << phiMin << "; phi max: " 
						<< phiMin + NphiBins*phiStep << "; phi step: " << phiStep << endl;
			out << "eta bins: " << NetaBins << "; eta min: " << etaMin << "; eta max: " 
						<< etaMin + NetaBins*etaStep << "; eta step: " << etaStep << endl;
			out << "invPT bins: " << NinvptBins << "; invPt min: " << invPtMin << "; invPt max: " 
						<< invPtMin + NinvptBins*invPtStep << "; invPt step: " << invPtStep << endl;
										
			//out << "ETA " << NetaBins << " " << etaMin << " " << etaStep << endl;
			//out << "INVPT " << NinvptBins << " " << invPtMin << " " << invPtStep << endl;
			
			//cout << "level: " << level << endl;
			
			if(level > 0) for(int i = 0; i != NphiBins; ++i) 
							for(int j = 0; j != NetaBins; ++j)
								for(int k = 0; k != NinvptBins; ++k){
								
									out << endl << "cell ijk: " << i << " " << j << " " << k << "      " 
										<< " phi:" << phiMin + i*phiStep << " "
										<< " eta:" << etaMin + j*etaStep << " "
										<< " invPt: " << invPtMin + k* invPtStep << endl;
										
									ArrElem[i][j][k].print(out);
								}
		}
		
		
		void write(ostream &out){
		
			out << TrainingPhase << endl;
			
			// parameters to check consistency of file with current HTM
			
			out << NphiBins << " " << phiMin << " " << phiStep << endl;
			out << NetaBins << " " << etaMin << " " << etaStep << endl;
			out << NinvptBins << " " << invPtMin << " " << invPtStep << endl;
			
				
			 for(int i = 0; i != NphiBins; ++i) 
				for(int j = 0; j != NetaBins; ++j)
					for(int k = 0; k != NinvptBins; ++k){										
						for(auto it = ArrElem[i][j][k].layerIndHitStat.begin(); it != ArrElem[i][j][k].layerIndHitStat.end(); ++it){  
								out << i << " " << j << " " << k << "  " << it -> first << " ";	// cell coordinates and layer  
								out << ArrElem[i][j][k].minLayers << " " << ArrElem[i][j][k].maxLayers << " ";
        						(it->second).write(out); // umin-umax
        						out << endl;
        				} 
					}
		}
		
		
		void clear(){
		
			for(int i = 0; i != NphiBins; ++i) 
				for(int j = 0; j != NetaBins; ++j)
					for(int k = 0; k != NinvptBins; ++k)												
						ArrElem[i][j][k].layerIndHitStat.clear();
        				 
					
		}
		
		
		int read(istream &in){
			
			int i,j,k,iLay;
			unsigned minLayers, maxLayers;
			double  angle, mean[3],low[3], high[3], u0L, u0H, u1L, u1H, u2L, u2H ;
			unsigned _NphiBins, _NetaBins, _NinvptBins;
			double _phiMin, _phiStep, _etaMin, _etaStep, _invPtMin, _invPtStep;
			int retcode;
						
			in >> retcode;
			
			// check consistency of file with current HTM
			
			double toll = 1.e-7;
			
			in >> _NphiBins >> _phiMin >> _phiStep;
			in >> _NetaBins >> _etaMin >> _etaStep;
			in >> _NinvptBins >> _invPtMin >>_invPtStep;
			
			if((_NphiBins-NphiBins) != 0) return -1;
			if(fabs(_phiMin-phiMin) > toll) return -2;
			if(fabs(_phiStep-phiStep) > toll) return -3;
			
			if((_NetaBins-NetaBins)!= 0) return -4;
			if(fabs(_etaMin-etaMin) > toll) return -5;
			if(fabs(_etaStep-etaStep) > toll) return -6;
			
			if((_NinvptBins-NinvptBins) != 0) return -7;
			if(fabs(_invPtMin-invPtMin)> toll) return -8;
			if(fabs(_invPtStep-invPtStep) > toll) return -9;
			
			// read file and fill HT array
		
			for(; in ;){
				in >> i >> j >> k >> iLay 
					>> minLayers >> maxLayers 
					>> mean[0] >> mean[1]>> mean[2] 
					>> low[0] >> low[1] >> low[2]  
					>> high[0] >> high[1] >> high[2]  
						>> angle >> u0L >> u1L >> u2L >> u0H >> u1H >> u2H ;
				//cout << "read "<< i << " " << j << " " << k << " " << iLay << " " 
				//	<< mean[0] << " " << mean[1] << " " << mean[2] << " " << angle << endl;
				if(i < 0 || i >= NphiBins) return -1;
				if(j < 0 || j >= NetaBins) return -1;
				if(k < 0 || k >= NinvptBins) return -1;
				Hitstat hs;
				hs.mean[0] = mean[0];
				hs.mean[1] = mean[1];
				hs.mean[2] = mean[2];
				hs.low[0] = low[0];
				hs.low[1] = low[1];
				hs.low[2] = low[2];
				hs.high[0] = high[0];
				hs.high[1] = high[1];
				hs.high[2] = high[2];
				hs.rotAngle = angle;
				hs.cosAngle = cos(angle);
				hs.sinAngle = sin(angle);
				hs.u0Min = u0L * 1.01;
				hs.u0Max = u0H * 1.01;
				hs.u1Min = u1L * 1.01;
				hs.u1Max = u1H * 1.01;
				hs.u2Min = u2L * 1.01;
				hs.u2Max = u2H * 1.01;
				
				ArrElem[i][j][k].layerIndHitStat[iLay]= hs;
				ArrElem[i][j][k].minLayers = minLayers;
				ArrElem[i][j][k].maxLayers = maxLayers;
			}
			if(retcode < 3) ++retcode;
			return retcode;
		
		} // end read
		
		
		// diagonalize each element of the array
		
		void diagonalize(){	
			
		 	for(int i = 0; i != NphiBins; ++i) 
				for(int j = 0; j != NetaBins; ++j)
					for(int k = 0; k != NinvptBins; ++k)	
									ArrElem[i][j][k].diagonalize();
		}
		
		// find cell i,j,k in HT array from track parameters
		
		int getCell(Track &t, int &i, int &j, int &k){
		
			 i = (t.phi - phiMin)/phiStep;
			 j = (t.eta - etaMin)/etaStep;
			 k = (t.invPt - invPtMin)/invPtStep;
			
			// check boundaries - return -1 if track parameters are out of bounds
			
			if(i < 0 || i >= NphiBins) return -1;
			if(j < 0 || j >= NetaBins) return -1;
			if(k < 0 || k >= NinvptBins) return -1;
			
			return 0;
		
		}
				
		
		// Train this Array with hit h coming from track t
		
		int train(Hit &h, Track &t){
		
			// find cell i,j,k in HT array from track parameters
			
			int i = (t.phi - phiMin)/phiStep;
			int j = (t.eta - etaMin)/etaStep;
			int k = (t.invPt - invPtMin)/invPtStep;
			
			// check boundaries - return -1 if track parameters are out of bounds
			
			if(i < 0 || i >= NphiBins) return -1;
			if(j < 0 || j >= NetaBins) return -1;
			if(k < 0 || k >= NinvptBins) return -1;
			
			// train the appropriate cell
			
			ArrElem[i][j][k].train(h);
			
			// fill max and min for num of layers traversed by this track
		     
        	unsigned nL = t.hitList.size();
        	if(nL < ArrElem[i][j][k].minLayers) ArrElem[i][j][k].minLayers = nL;
        	if(nL > ArrElem[i][j][k].maxLayers) ArrElem[i][j][k].maxLayers = nL;
   		
				// Fill the 3D plot just for the one special test cell
				
				if(i==plotBinX && j==plotBinY && k==plotBinZ){
					int lay = h.iLayer;
					if (h.hitType == 'D') lay += dg.nBarrels;
					if(lay == plotLay) {
						// This is the right cell we want to plot
						if(H3D_HTcellx->GetEntries()==0){
							// Set histogram boundaries
							double xLow = ArrElem[i][j][k].layerIndHitStat[lay].low[0];
							double xHigh = ArrElem[i][j][k].layerIndHitStat[lay].high[0];
							double yLow = ArrElem[i][j][k].layerIndHitStat[lay].low[1];
							double yHigh = ArrElem[i][j][k].layerIndHitStat[lay].high[1];
							double zLow = ArrElem[i][j][k].layerIndHitStat[lay].low[2];
							double zHigh = ArrElem[i][j][k].layerIndHitStat[lay].high[2];
							H3D_HTcellx->GetXaxis()->SetLimits(xLow,xHigh);
							H3D_HTcellx->GetYaxis()->SetLimits(yLow,yHigh);
							H3D_HTcellx->GetZaxis()->SetLimits(zLow,zHigh);	
							
							xLow = ArrElem[i][j][k].layerIndHitStat[lay].u0Min;
							xHigh = ArrElem[i][j][k].layerIndHitStat[lay].u0Max;
							yLow = ArrElem[i][j][k].layerIndHitStat[lay].u1Min;
							yHigh = ArrElem[i][j][k].layerIndHitStat[lay].u1Max;
							zLow = ArrElem[i][j][k].layerIndHitStat[lay].u2Min;
							zHigh = ArrElem[i][j][k].layerIndHitStat[lay].u2Max;
							H3D_HTcellu->GetXaxis()->SetLimits(xLow,xHigh);
							H3D_HTcellu->GetYaxis()->SetLimits(yLow,yHigh);
							H3D_HTcellu->GetZaxis()->SetLimits(zLow,zHigh);	
							
							//cout << xLow << " " << xHigh << " " << yLow << " " << yHigh
							//<< " " << zLow << " " << zHigh << endl;
						}
						H3D_HTcellx->Fill(h.x1,h.x2,h.t);
						H3D_HTcellu->Fill(h.u0,h.u1,h.u2);
					}
				}
		
			
			return 0;
			
		} // end train
		
		
		void reset(){
		
			for(unsigned iPhi = 0; iPhi != NphiBins; ++iPhi)
				for(unsigned iEta = 0; iEta != NetaBins; ++iEta)
					for(unsigned iInvpt = 0; iInvpt != NinvptBins; ++iInvpt)
						ArrElem[iPhi][iEta][iInvpt].reset();
			
			return;
		
		
		}
		
		
		
		int fill(Hit &h, int mode = 0){
		
			int iPhi1 = 0;
			int iPhi2 = NphiBins;
			int iEta1 = 0;
			int iEta2 = NetaBins -1;
			int iInvpt1 = 0;
			int iInvpt2 = NinvptBins;
		
			
			if(mode == 1){
			
				double r,z;
				if(h.hitType == 'B'){
					z = h.x2;
					r = dg.B[h.iLayer].r;	
				}
				else{
					z = dg.D[h.iLayer].z;
					r = h.x2;	
				}		
				double eta = asinh(z/r);

				// index of the hit eta bin
				int iEtaMap = (eta-etaMin)/etaStep;

				int errEta = 2;
				iEta1 = iEtaMap - errEta;
				iEta2 = iEtaMap + errEta +1;
				
				if(iEta2 < 0) return 0;
				if(iEta1 < 0) iEta1 = 0;
				if(iEta1 > NetaBins-1) iEta1 = NetaBins - 1;
				if(iEta2 > NetaBins) iEta2 = NetaBins;
				
			} // end mode = 1
			
				for(unsigned iPhi = iPhi1; iPhi != iPhi2; ++iPhi){
						for(unsigned iEta = iEta1; iEta != iEta2; ++iEta){
							for(unsigned iInvpt = iInvpt1; iInvpt != iInvpt2; ++iInvpt){	
								int strike = ArrElem[iPhi][iEta][iInvpt].fill(h, mode);			
								if(strike == 0) {
									Hitstat* thisStat = &(ArrElem[iPhi][iEta][iInvpt].layerIndHitStat[h.layerInd]);
									thisStat->hitLayer = true;
									thisStat->hitIDList.push_back(h.ID);
								}							
							} // end loop on invpt	
						} // end loop on eta
				} // end loop on phi	
			
			
			return 0;
			
		}// end fill
			
			
		unsigned getBestCell(unsigned &iPhi_, unsigned &iEta_, unsigned &iInvpt_){

			unsigned maxHits = 0;
			for(unsigned iPhi = 0; iPhi != NphiBins; ++iPhi)
					for(unsigned iEta = 0; iEta != NetaBins; ++iEta)
						for(unsigned iInvpt = 0; iInvpt != NinvptBins; ++iInvpt){
							unsigned nHits = ArrElem[iPhi][iEta][iInvpt].nHitLayers;

							if(nHits > maxHits){
								maxHits = nHits;
								iPhi_ = iPhi; iEta_ = iEta; iInvpt_ = iInvpt;
							}

						}

			return maxHits;
		    			
		    } // end getBestCell
		
		struct Pars{unsigned iPhi; unsigned iEta; unsigned iInvpt; unsigned nLayers;};
		vector<Pars> cellCandidateList;
		
			
		unsigned getCellCandidates(){
		
			cellCandidateList.clear();
			for(unsigned iPhi = 0; iPhi != NphiBins; ++iPhi)
				for(unsigned iEta = 0; iEta != NetaBins; ++iEta)
					for(unsigned iInvpt = 0; iInvpt != NinvptBins; ++iInvpt){
						unsigned nLayers = ArrElem[iPhi][iEta][iInvpt].nHitLayers;
						if(nLayers >= ArrElem[iPhi][iEta][iInvpt].minLayers){
							Pars p;
							p.iPhi = iPhi;
							p.iEta = iEta;
							p.iInvpt = iInvpt;
							p.nLayers = nLayers;
							cellCandidateList.push_back(p);					
						}
					}
		    return (unsigned)cellCandidateList.size();				
		    			
		} // end getCellCandidates
		
		
		void printCellCandidateList(ostream &out){
		
			for(unsigned i = 0; i != cellCandidateList.size(); ++i){
				unsigned iPhi = cellCandidateList[i].iPhi;
				unsigned iEta = cellCandidateList[i].iEta;
				unsigned iInvpt = cellCandidateList[i].iInvpt;
				unsigned nLayers = cellCandidateList[i].nLayers;	
				out << "Candidate " << i << ": " << iPhi << " " << iEta << " " << iInvpt << " nLayers: " << nLayers << endl;
				HTArrayElement elem = ArrElem[iPhi][iEta][iInvpt];
				elem.printHits(out);
			}
		
		
		} // end printCellCandidateList
    
    void writeCellCandidateList(ostream &out, vector<Hit> & hitList){
			out << cellCandidateList.size() << endl;
			for(unsigned i = 0; i != cellCandidateList.size(); ++i){
				unsigned iPhi = cellCandidateList[i].iPhi;
				unsigned iEta = cellCandidateList[i].iEta;
				unsigned iInvpt = cellCandidateList[i].iInvpt;
				unsigned nLayers = cellCandidateList[i].nLayers;	
				out << iPhi << " " << iEta << " " << iInvpt << " " << nLayers << endl;
				HTArrayElement elem = ArrElem[iPhi][iEta][iInvpt];
				elem.writeHits(out, hitList);
			}
		
		
		} // end writeCellCandidateList
    
    }; // end HTArray
         
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
