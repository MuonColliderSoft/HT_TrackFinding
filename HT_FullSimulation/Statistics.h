//
//  Statistics.h
//  LinearizedTrackFitting
//
//  Created by Luciano Ristori on 5/5/14.
//  Copyright (c) 2014 Luciano Ristori. All rights reserved.
//

#ifndef __LinearizedTrackFitting__Statistics__
#define __LinearizedTrackFitting__Statistics__

#include <iostream>

class Statistics {
    
public:
    
    long int nEntries;
    double mean;
    double variance;
    double max;
    double min;
    
    Statistics();
    
    void Fill(double x);
    
    long int GetEntries();
    double GetMean();
    double GetVariance();
    double GetSigma();
    double GetMax();
    double GetMin();
};

#endif /* defined(__LinearizedTrackFitting__Statistics__) */
