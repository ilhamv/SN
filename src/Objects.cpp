#include "Objects.h"

#include <algorithm>
#include <iostream>


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

void Region::reset_tau()
{
    r_tau = r_dz * M->SigmaT();
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
