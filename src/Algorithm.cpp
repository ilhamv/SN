#include <cmath>

#include "Algorithm.h"

//==========================================================================
// Miscellany
//==========================================================================

double norm_2( const std::vector<double> v )
{
    double norm = 0.0;
    for( int i = 0; i < v.size(); i++ ){
        norm += v[i]*v[i];
    }
    norm = std::sqrt(norm);
    return norm;
}

