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
void Material::set_tilde( const double zeta )
{
    m_SigmaS += zeta * m_nuSigmaF;
    m_SigmaA = m_SigmaT - m_SigmaS;
}
void Material::revert_tilde( const double zeta )
{
    m_SigmaS -= zeta * m_nuSigmaF;
    m_SigmaA = m_SigmaT - m_SigmaS;
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
void BCVacuum::set_boundary( std::vector<double>& psi, const int a, 
                             const int b )
{
    std::fill(psi.begin()+a, psi.begin()+b, 0.0);
}

// Reflective: Reverse psi
void BCReflective::set_boundary( std::vector<double>& psi, const int a, 
                                 const int b )
{
    for( int i = a; i < b; i++ ){
        psi[i] = psi[psi.size()-1-i];
    }
}

// Isotropic
void BCIsotropic::set_boundary( std::vector<double>& psi, const int a, 
                                const int b )
{
    std::fill(psi.begin()+a, psi.begin()+b, magnitude);
}

// Mono Directional
void BCMonoDirectional::set_boundary( std::vector<double>& psi, const int a, 
                                      const int b )
{
    std::fill(psi.begin()+a, psi.begin()+b, 0.0);
    psi[idx] = val;
}

// Linear
void BCLinear::set_boundary( std::vector<double>& psi, const int a, 
                                      const int b )
{
    int idx = 0;
    for( int i = a; i < b; i++ ){
        psi[i] = psi_b[idx]; idx++;
    }
}
