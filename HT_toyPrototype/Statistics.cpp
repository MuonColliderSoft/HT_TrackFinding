//
//  Statistics.cpp
//  LinearizedTrackFitting
//
//  Created by Luciano Ristori on 5/5/14.
//  Copyright (c) 2014 Luciano Ristori. All rights reserved.
//

#include <cmath>
#include "Statistics.h"


Statistics::Statistics() {
    
    nEntries = 0;;
    mean = 0.;
    variance = 0.;
    max = 0.;
    min = 0.;
    
};

void Statistics::Fill(double x){
    ++ nEntries;
    mean += (x - mean)/nEntries;
    if(nEntries >=2) variance += (x - mean)*(x - mean)/(nEntries-1) - variance/nEntries;
    if(nEntries == 0) {
        max = x;
        min = x;
    }
    else if(x > max) max = x;
    else if(x < min) min = x;
}

long int Statistics::GetEntries(){
    return nEntries;
};
double Statistics::GetMean() {
    return mean;
};
double Statistics::GetVariance(){
    return variance;
};
double Statistics::GetSigma(){
    return sqrt(variance);
    
};

double Statistics::GetMax(){
    return max;
};

double Statistics::GetMin(){
    return min;
};


