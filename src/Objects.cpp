#include "Objects.h"

#include <algorithm>
#include <iostream>
#include <cmath>


//==============================================================================
// Material
//==============================================================================

void Material::time_augment( const double aug )
{
    m_SigmaT += aug;
    m_SigmaA += aug;
}
void Material::revert_augment( const double aug )
{
    m_SigmaT -= aug;
    m_SigmaA -= aug;
}

//==============================================================================
// Region
//==============================================================================

void Region::reset()
{
    r_tau = r_dz * M->SigmaT();
}
        
void Region::set_alpha( const std::vector<double>& mu, const std::string type )
{
    const double N = mu.size();
    r_alpha.resize(N);

    for( int n = 0; n < N; n++ ){
        if( type == "DD" ) { r_alpha[n] = 0.0; } 
        else if ( type == "SC" ){ 
            r_alpha[n] = 1.0 / ( std::tanh( 0.5 * r_tau / mu[n] ) )
                         - ( 2.0 * mu[n] / r_tau );
        }
    }
}

//==============================================================================
// Boundary Condition
//==============================================================================

// Vacuum: Zero out psi
void BCVacuum::set_boundary( std::vector<double>& psi )
{
    std::fill(psi.begin(), psi.end(), 0.0);
}

// Reflective: Reverse psi
void BCReflective::set_boundary( std::vector<double>& psi )
{
    std::reverse(psi.begin(), psi.end());
}

// Isotropic
void BCIsotropic::set_boundary( std::vector<double>& psi )
{
    std::fill(psi.begin(), psi.end(), magnitude);
}

// Mono Directional
void BCMonoDirectional::set_boundary( std::vector<double>& psi )
{
    std::fill(psi.begin(), psi.end(), 0.0);
    psi[idx] = val;
}
