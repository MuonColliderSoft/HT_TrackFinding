//
//  Gauss.h
//

double gauss (std::mt19937 &generator){
    
    double sum = 0.;
    for(int i = 0; i !=11; ++i) sum += distribution(generator);
    return sum - 6.;
};

