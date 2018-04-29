#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>

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

double norm_2( const std::vector<double> v, 
               const std::vector<std::shared_ptr<Region>>& mesh )
{
    double norm = 0.0;
    for( int i = 0; i < v.size(); i++ ){
        norm += v[i]*v[i] * mesh[i]->dz();
    }
    norm = std::sqrt(norm);
    return norm;
}

double wielandt_shift( const double ws_scale, const double ws_subtract,
                       const double ws_min, const double lambda_old,
                       const double lambda_new )
{
    return std::max( ws_scale * lambda_new - 
                     ws_subtract * std::abs( lambda_new - lambda_old ),
                     ws_min );
}
