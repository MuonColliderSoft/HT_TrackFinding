//
//  HTArray.cpp
//  MuonColliderToy
//
//  Created by Luciano Ristori on 5/23/22
//


 extern int TrainingPhase;
 extern bool Special;
 extern bool verbose;

// global variables	
 	
vector<Hit> fitHitList; // list of all hits to be used in one fit
vector< vector <Hit> > combinations; // all hit combinations to be fitted in one candidate


double xi_to_phi(double xi){return par.HTA_t_phi - par.HTA_t_deltaPhi + (2*par.HTA_t_deltaPhi/(double)HTA_NphiBins)*xi;}
double xj_to_eta(double xj){return par.HTA_t_eta - par.HTA_t_deltaEta + 2*(par.HTA_t_deltaEta/(double)HTA_NetaBins)*xj;}
double xk_to_invPt(double xk){return par.HTA_t_invPt_min + (par.HTA_t_invPt_max - par.HTA_t_invPt_min)/(double)HTA_NinvptBins*xk;}

double phi_to_xi(double phi){return (phi - par.HTA_t_phi + par.HTA_t_deltaPhi)/par.HTA_t_deltaPhi*(double)HTA_NphiBins/2.;}
double eta_to_xj(double eta){return (eta - par.HTA_t_eta + par.HTA_t_deltaEta)/par.HTA_t_deltaEta*(double)HTA_NetaBins/2.;}
double invPt_to_xk(double invPt){return invPt/(par.HTA_t_invPt_max - par.HTA_t_invPt_min)*(double)HTA_NinvptBins;}




double xxxxchi2Func(const double *x){
	 	   	       	
	 double invPt = xk_to_invPt(x[0]);
	 double eta = xj_to_eta(x[1]);
	 double phi = xi_to_phi(x[2]);
	 double z0 = x[3];
	 double t0 = x[4];
	 
	int particleType = 2; //  pion

	Track tFit(particleType, 0., 0., z0, t0, invPt, eta, phi); // specific track constructor

	double chi2 = 0.;
	for(int iHit = 0; iHit != fitHitList.size(); ++iHit){
		chi2 += tFit.hitChi2(dg, fitHitList[iHit]);
	}

	return chi2;

}// end chi2Func


double chi2MassFunc(const double *x){
	 	   	       	
	 double invPt = xk_to_invPt(x[0]);
	 double eta = xj_to_eta(x[1]);
	 double phi = xi_to_phi(x[2]);
	 double z0 = x[3];
	 double t0 = x[4];
	 double mass = x[5];
	 

	Track tFit(mass, 0., 0., z0, t0, invPt, eta, phi); // specific track constructor

	double chi2 = 0.;
	for(int iHit = 0; iHit != fitHitList.size(); ++iHit){
		chi2 += tFit.hitChi2(dg, fitHitList[iHit]);
	}

	return chi2;

}// end chi2Func



 	

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

    class Hitstat {
    
    	// This is the class that contains all that is needed for each layer in each HT cell
    
    	    
		public:
		
		bool hitLayer; // this layer has been hit in this HT element	
		vector<unsigned> hitIDList; // list of hit ID's (indices in hit vector in event)
    	unsigned nHits;
    	
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
    		nHits = 0;
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
    		nHits = 0;
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
       	
		int i_index;
		int j_index;   
		int k_index; 
    		
		double imeanPhi;
		double jmeanEta;
		double kmeanInvPt;
		
		bool thisCellDone = false;
               
        unsigned nBarrels = dg.nBarrels;       
        unsigned nHitLayers	; // number of layers hit for current event in this element
             
        // list of layers for this element implemented as a map
        map<int, Hitstat> layerIndHitStat;
        map<int, Hitstat>::iterator it;
        
        // next two params calculated in training phase and loaded at initialization
        unsigned minLayers = 10000; // min number of layers hit by a track from this param space
        unsigned maxLayers = 0;// max number of layers hit by a track from this param space
        
       // train this array element
        // during TrainigPhase 1, add this layer to the list for this element if missing
        
        void train(Hit &h){
        
        	if(TrainingPhase == 1)layerIndHitStat[h.layerInd].train(h);
        	else if(layerIndHitStat.find(h.layerInd) != layerIndHitStat.end())
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
         
         	thisCellDone = false;
         
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
			
				/*if(i_index == 29 && j_index == 97 && k_index == 9) {
					Special = true;			
					cout << "cell [" << i_index <<  "][" << j_index << "][" << k_index << "]" << endl;
					cout << "hit "; h.print(cout);			
				}*/
									
				double u0 =  h.x1 - stat.mean[0];
				if(Special) cout << u0 << " " << stat.u0Min << " " 	<< stat.u0Max << endl;	
				if(u0 < stat.u0Min || u0 > stat.u0Max) {Special = false; return -1;}
					if(Special)cout << "good u0" << endl;
				double u1 = (h.x2 - stat.mean[1]) * stat.cosAngle - (h.t - stat.mean[2]) * stat.sinAngle;
				if(Special) cout << u1 << " " << stat.u1Min << " " 	<< stat.u2Max << endl;		
				if(u1 < stat.u1Min || u1 > stat.u1Max) {Special = false; return -2; }
					if(Special)cout << "good u1"<< endl;
				double u2 = (h.x2 - stat.mean[1]) * stat.sinAngle + (h.t - stat.mean[2]) * stat.cosAngle;	
				if(Special) cout << u2 << " " << stat.u2Min << " " 	<< stat.u2Max << endl;		
				if(u2 < stat.u2Min || u2 > stat.u2Max) { Special = false; return -3;}
					if(Special)cout << "good u2"<< endl;
						
				if(!stat.hitLayer)++nHitLayers;
				++(it->second.nHits);
				Special = false;
				return 0;
					
			} 
			
			else return -4; 
			  
        }// end fill
        
        
        void printCandidate(Event &event, std::ostream &out){
        
        	// iterate on all layers 
			  
			for(it = layerIndHitStat.begin(); it != layerIndHitStat.end(); ++it){ 
			
					// iterate on all hits in this layer
				 
					int nHits = (it->second).hitIDList.size();
					for(int ih = 0; ih != nHits; ++ih){
						Hit h = event.hitList[(it->second).hitIDList[ih]];
						cout << "L" << it->first << " ";
						h.print(cout);
												
					}
								
			} 
        
        
        } // end printCandidate
        
        
        
        int fitCandidate(Event &event, double &chi2, double &phi, double &eta, double &invPt, double &z0, double &t0, double &mass, unsigned &nLayers){
        
        	struct Result {
        		double chi2;
        		double phi;
        		double eta;
        		double invPt; 
        		double z0;
        		double t0;
        		double mass; 
        		unsigned nLayers;
        	};
        	
        	Result bestFitRes;
        	int bestFitComb;
        
        	// create list of combination of hits to be fitted
        	
			vector< vector <Hit> > combinations2; // all hit combinations to be fitted in one candidate
			vector< vector <Hit> > *comb1 = &combinations;
			vector< vector <Hit> > *comb2 = &combinations2;
			vector< vector <Hit> > *temp;
			vector<Hit> hitList;
			
			// create one single empty combination
			
			hitList.clear();
			comb1->clear();
			comb2->clear();
			comb1->push_back(hitList);
			
			// iterate on all layers 
			  
			for(it = layerIndHitStat.begin(); it != layerIndHitStat.end(); ++it){ 
			
				int nHits = (it->second).hitIDList.size();
				if(nHits == 0) continue;// skip possible empty layers
				
				// iterate on all previous combinations
				int nComb = comb1->size();
				for(int iC = 0; iC != nComb; ++iC){ 
				
					// iterate on all hits in this layer
				 
					for(int ih = 0; ih != nHits; ++ih){
						Hit h = event.hitList[(it->second).hitIDList[ih]];
							hitList = (*comb1)[iC];
							hitList.push_back(h);
							comb2->push_back(hitList);
												
					}
				(*comb1)[iC].clear();						
				}
			
			comb1->clear();
			temp = comb1;
			comb1 = comb2;
			comb2 = temp;
						
			} 
			
			combinations = *comb1;
			
			//nHits = combinations[0].size();
			int nCombs = combinations.size();
			
			if(verbose) cout << nCombs << " combinations" << endl;
        	
		
			// Minimize with Root Minimizer
	
			// create minimizer giving a name and a name (optionally) for the specific
			// algorithm
			// possible choices are: 
			//     minName                  algoName
			// Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
			//  Minuit2                     Fumili2
			//  Fumili
			//  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS, 
			//                              BFGS2, SteepestDescent
			//  GSLMultiFit
			//   GSLSimAn xxx problems
			//   Genetic
		
			const char * minName = "Minuit2";
			const char *algoName = "";
		
			ROOT::Math::Minimizer* min = 
			  ROOT::Math::Factory::CreateMinimizer(minName, algoName);
			
			// loop on all combinations for this candidate
			
			bool oneGoodFit = false; // did we find at least one good hit combination?
			bool firstGoodFit = true; // is this the first good fit for min chi2?
			
			for(int iComb = 0; iComb != nCombs; ++iComb){ 
			
				fitHitList = combinations[iComb];// pick the correct hit combination for chi2Funcx
				
				if(verbose) {
					cout << "combination " << iComb;			
					cout << " - ";				
					for(unsigned iH = 0; iH != fitHitList.size(); ++iH){
						cout << fitHitList[iH].trackInd << " ";
					}
					cout << endl;
					for(unsigned iH = 0; iH != fitHitList.size(); ++iH){
						fitHitList[iH].print(cout);
					}
				}
				
				// set tolerance , etc...
				min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
				min->SetMaxIterations(1000);  // for GSL 
				min->SetTolerance(0.001);
				min->SetPrintLevel(0);

				// create function wrapper for minmizer
				// a IMultiGenFunction type 
				ROOT::Math::Functor f(&chi2MassFunc,6); 
				double step[6] = {0.1,0.1,0.1,0.1,0.1,0.001};
				double variable[6] = {kmeanInvPt,jmeanEta,imeanPhi,0.0,0.0,0.139570};
		   		
		   		
				 min->SetFunction(f);
				 
				 
				
				 
				// Set the free variables to be minimized!
		
				min->SetVariable(0,"invPt",variable[0], step[0]);
				min->SetVariable(1,"Eta",variable[1], step[1]);
				min->SetVariable(2,"Phi",variable[2], step[2]);
				min->SetVariable(3,"z0",variable[3], step[3]);
				min->SetVariable(4,"t0",variable[4], step[4]);
				min->SetVariable(5,"mass",variable[5], step[5]);
				
				 // do the first minimization with  mass fixed
				
				min->FixVariable(5);			
		
				min->Minimize(); 
		
				if(par.gen_massFit){
				
					// do the second minimization with free mass and t0
				
					min->ReleaseVariable(5);
					min->FixVariable(0);
					min->FixVariable(1);
					min->FixVariable(2);
					min->FixVariable(3);
							
				min->Minimize(); 
				}
		
				const double *res = min->X();
						
				// return parameters
				
				chi2 = min->MinValue(); 
				phi = xi_to_phi(res[2]);
				eta = xj_to_eta(res[1]);
				invPt = xk_to_invPt(res[0]);
				z0 = res[3];				
				t0 = res[4];
				mass = res[5];
				 
				
				if(verbose) cout << "chi2: " << chi2 << endl;

				
				oneGoodFit = true;
				
				if(firstGoodFit || bestFitRes.chi2 > chi2){
					firstGoodFit = false;
					bestFitRes.chi2 = chi2;
					bestFitRes.phi = phi;
					bestFitRes.eta = eta;
					bestFitRes.invPt = invPt;
					bestFitRes.z0 = z0;
					bestFitRes.t0 = t0;
					bestFitRes.mass = mass;
					bestFitComb = iComb;				
									
				}
			
			}// end loop on combinations
			
			if(oneGoodFit){
				nLayers = combinations[bestFitComb].size();
		
				chi2 = bestFitRes.chi2;
				phi = bestFitRes.phi;
				eta =bestFitRes.eta;
				invPt = bestFitRes.invPt;
				z0 = bestFitRes.z0;
				t0 = bestFitRes.t0;
				mass = bestFitRes.mass;
				
				if(verbose) cout << "bestFitComb: " << bestFitComb << " chi2: " << chi2 << endl;
				
				return 0;
			}
			else return -1;
	
        }// end fitCandidate
        
    }; // end HTArrayElement
    
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
    
    class HTArray {
               
    // This is the whole Hough Transform Array with parameters and accessories
  
    	public:
    	
    	//ofstream outfile;	
    	
    	// some histograms of inside quantities
	
		TH1D *HDeltaEta;
		TH1D *HDeltaPhi;
		TH1D *HDeltaPhi2;
		TH1D *HInvptTimesR; 
		TH2I *HDeltaPhiVsInvptTimesR;
		
		
		// histograms of a single HT cell content as an example
		
		TH3I *H3D_HTcellx;
		TH3I *H3D_HTcellu;
		TH1D *Hcellx1;
		TH1D *Hcellx2;
		TH1D *Hcellt;
		TH1D *Hcellu0;
		TH1D *Hcellu1;
		TH1D *Hcellu2;

    
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
                
        double signedCurvaturePerInvPt;
        
        double getCell_Phi(int i){
			return phiMin + ((double)i+0.5)*phiStep;
		}
			
		double getCell_Eta(int j){
			return etaMin + ((double)j+0.5)*etaStep;
		}
		
		double getCell_InvPt(int k){
			return invPtMin + ((double)k+0.5)*invPtStep;
		}	
				
         
        double get_Phi(double xi){
			return phiMin + xi*phiStep;
		}
		
		double get_Eta(double xj){
			return etaMin + xj*etaStep;
		}
		
		double get_InvPt(double xk){
			return invPtMin + xk*invPtStep;
		}	
		
		
		HTArray(){ // constructor
		
			phiMin = par.HTA_t_phi - par.HTA_t_deltaPhi;
			phiStep = 2*par.HTA_t_deltaPhi/double(NphiBins);
			etaMin = par.HTA_t_eta - par.HTA_t_deltaEta;
			etaStep = 2*par.HTA_t_deltaEta/double(NetaBins);
			invPtMin = par.HTA_t_invPt_min;
			invPtStep = (par.HTA_t_invPt_max - par.HTA_t_invPt_min)/double(NinvptBins);
			
			// constant for phi cell prediction
			
			signedCurvaturePerInvPt = 3.e-4*par.magneticField;
		
						
			 for(int i = 0; i != NphiBins; ++i) 
				for(int j = 0; j != NetaBins; ++j)
					for(int k = 0; k != NinvptBins; ++k){										
						ArrElem[i][j][k].imeanPhi = i + 0.5;										
						ArrElem[i][j][k].jmeanEta = j + 0.5;											
						ArrElem[i][j][k].kmeanInvPt = k + 0.5;											
						ArrElem[i][j][k].i_index = i;										
						ArrElem[i][j][k].j_index = j;											
						ArrElem[i][j][k].k_index = k;	
					}		
		
		}
		
		void initHists(){
			
			HDeltaEta = new TH1D("DeltaEta","DeltaEta",21,-10.5,+10.5);
			HDeltaPhi = new TH1D("DeltaPhi","DeltaPhi",51,-10.5,+40.5);
			HDeltaPhi2 = new TH1D("DeltaPhi2","DeltaPhi2",21,-10.5,+10.5);
			HInvptTimesR = new TH1D("InvptTimesR","InvptTimesR",1000, 0.,0.3);
			HDeltaPhiVsInvptTimesR = new TH2I("DeltaPhiVsInvptTimesR","DeltaPhiVsInvptTimesR",100,0.,0.3,510,-10.5,+40.5);					 
			H3D_HTcellx = new TH3I("H3D_HTcellx","H3D_HTcellx; x1; x2; t",40,0.,1.,40,0.,1.,40,0.,1.);
			H3D_HTcellu = new TH3I("H3D_HTcellu","H3D_HTcellu; x1; u1; u2",40,0.,1.,40,0.,1.,40,0.,1.);
			Hcellx1 = new TH1D("Hcellx1","Hcellx1",100,0.,1.);
			Hcellx2 = new TH1D("Hcellx2","Hcellx2",100,0.,1.);
			Hcellt = new TH1D("Hcellt","Hcellt",100,0.,1.);
			Hcellu0 = new TH1D("Hcellu0","Hcellu0",100,0.,1.);
			Hcellu1 = new TH1D("Hcellu1","Hcellu1",100,0.,1.);
			Hcellu2 = new TH1D("Hcellu2","Hcellu2",100,0.,1.);
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
		
		int getCell_i(double xPhi, int &i){
			i = (xPhi - phiMin)/phiStep;
			if(i < 0 || i >= NphiBins) return -1;
			else return 0;	
		}
		
		int getCell_j(double xEta, int &j){
			j = (xEta - etaMin)/etaStep;
			if(j < 0 || j >= NetaBins) return -1;
			else return 0;	
		}
		
		int getCell_k(double xInvPt, int &k){
			k = (xInvPt - invPtMin)/invPtStep;
			if(k < 0 || k >= NinvptBins) return -1;
			else return 0;	
		}
		
	
		
		// Train this Array with hit h coming from track t
		
		int train(Hit &h, Track t, unsigned i, unsigned j, unsigned k){
				
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
						// This is the special cell we want to plot
						
						if(H3D_HTcellx->GetEntries()==0){// init hist limits (first fill only)
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
							
							Hcellx1->GetXaxis()->SetLimits(xLow,xHigh);
							Hcellx2->GetXaxis()->SetLimits(yLow,yHigh);
							Hcellt->GetXaxis()->SetLimits(zLow,zHigh);
							
							xLow = ArrElem[i][j][k].layerIndHitStat[lay].u0Min;
							xHigh = ArrElem[i][j][k].layerIndHitStat[lay].u0Max;
							yLow = ArrElem[i][j][k].layerIndHitStat[lay].u1Min;
							yHigh = ArrElem[i][j][k].layerIndHitStat[lay].u1Max;
							zLow = ArrElem[i][j][k].layerIndHitStat[lay].u2Min;
							zHigh = ArrElem[i][j][k].layerIndHitStat[lay].u2Max;
							H3D_HTcellu->GetXaxis()->SetLimits(xLow,xHigh);
							H3D_HTcellu->GetYaxis()->SetLimits(yLow,yHigh);
							H3D_HTcellu->GetZaxis()->SetLimits(zLow,zHigh);	
							
							Hcellu0->GetXaxis()->SetLimits(xLow,xHigh);
							Hcellu1->GetXaxis()->SetLimits(yLow,yHigh);
							Hcellu2->GetXaxis()->SetLimits(zLow,zHigh);
						
						}// end init hist limits
						
						H3D_HTcellx->Fill(h.x1,h.x2,h.t);
						H3D_HTcellu->Fill(h.u0,h.u1,h.u2);					
						Hcellx1->Fill(h.x1);
						Hcellx2->Fill(h.x2);
						Hcellt->Fill(h.t);
						Hcellu0->Fill(h.u0);
						Hcellu1->Fill(h.u1);
						Hcellu2->Fill(h.u2);
						
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
			int iEta2 = NetaBins;
			int iInvpt1 = 0;
			int iInvpt2 = NinvptBins;
		
			
			
			
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

				// indexes of the hit coordinates
						
				int iPhiMap;
				int iPhiMap2;
				
				
				int iEtaMap; 
				getCell_j(eta, iEtaMap);
				

			if(mode == 1 || mode == 2){
	
				int errEta = par.gen_errEta;
				iEta1 = iEtaMap - errEta;
				iEta2 = iEtaMap + errEta + 1;
				
				if(iEta2 < 1) return 0;
				if(iEta1 < 0) iEta1 = 0;
				if(iEta1 > NetaBins-1) iEta1 = NetaBins - 1;
				if(iEta2 > NetaBins) iEta2 = NetaBins;
				
			} // end mode = 1 or 2
			
				//h.print(cout);
		
				for(unsigned iEta = iEta1; iEta != iEta2; ++iEta){
					for(unsigned iInvpt = iInvpt1; iInvpt != iInvpt2; ++iInvpt){
					
						//cout << iEta << " "<< iInvpt << endl;
					
						getCell_i(h.realPhi(dg), iPhiMap);
						getCell_i(h.realPhi(dg) - 0.5*signedCurvaturePerInvPt*getCell_InvPt(iInvpt)*r, iPhiMap2);
						
						if(mode == 2){
							int errPhi = par.gen_errPhi;
							//cout << "errPhi = " << errPhi << endl;
							iPhi1 = iPhiMap2 - errPhi;
							iPhi2 = iPhiMap2 + errPhi +1;
				
							if(iPhi2 < 1) continue;
							if(iPhi1 < 0) iPhi1 = 0;
							if(iPhi1 > NphiBins-1) iPhi1 = NphiBins - 1;
							if(iPhi2 > NphiBins) iPhi2 = NphiBins;
							
							//cout << iPhi1 << " " << iPhi2 << endl;
				
						} // end mode = 2
						
						
					
						for(unsigned iPhi = iPhi1; iPhi != iPhi2; ++iPhi){
						//for(unsigned iPhi = 0; iPhi != NphiBins; ++iPhi){
						
							int strike = ArrElem[iPhi][iEta][iInvpt].fill(h, mode);			
							if(strike == 0) { //cout << "strike iPhi = " << iPhi << endl;
								Hitstat* thisStat = &(ArrElem[iPhi][iEta][iInvpt].layerIndHitStat[h.layerInd]);
								thisStat->hitLayer = true;
								thisStat->hitIDList.push_back(h.ID);
								
								// insert here code to build candidate list on the fly ////////////
							
								int deltaEta = iEtaMap - iEta;
								HDeltaEta->Fill(deltaEta);
							
								int deltaPhi = iPhiMap - iPhi;
								HDeltaPhi->Fill(deltaPhi);
								
								double xCoord = 0.5*signedCurvaturePerInvPt*getCell_InvPt(iInvpt)*r;
								HInvptTimesR->Fill(xCoord);
								HDeltaPhiVsInvptTimesR->Fill(xCoord,deltaPhi);
															
								int deltaPhi2 = iPhiMap2 - iPhi;			
								HDeltaPhi2->Fill(deltaPhi2);
								
								//cerr << iPhi << " " << iPhiMap << " " <<  iPhiMap2 << endl; // ********
							
									
								}							
							} // end loop on phi	
						} // end loop on pt
				} // end loop on eta	
			
			
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
						//if(nLayers >= (ArrElem[iPhi][iEta][iInvpt].minLayers))
						if(nLayers >= (ArrElem[iPhi][iEta][iInvpt].minLayers-1))					
							if(nLayers >= par.gen_minLayersForFit) {
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
