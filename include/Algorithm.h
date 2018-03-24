#ifndef _ALGORITHM_H
#define _ALGORITHM_H

#include <vector>
#include <sstream>
#include <memory>
#include <iterator>


//==========================================================================
// Parser tools
//==========================================================================

// Parsing space(s) delimited strings into vector
template<class T>
std::vector<T> parse_vector(std::string const& pointLine)
{
    std::istringstream iss(pointLine);

    return std::vector<T>{ std::istream_iterator<T>(iss), 
                           std::istream_iterator<T>()
                         };
}

// Find object pointer by id number
template< typename T >
std::shared_ptr<T> find_by_id( const std::vector<std::shared_ptr<T>>& vec,
                               const int id )
{
    for ( auto& v : vec ){
	if ( v->id() == id ) { return v; }
    }
    return nullptr;
}


//==========================================================================
// Miscellany
//==========================================================================

double norm_2( const std::vector<double> v );
void legendre_compute_glr ( int n, double x[], double w[] );

#endif // _ALGORITHM_H
