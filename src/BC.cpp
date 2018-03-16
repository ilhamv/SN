#include "BC.h"

#include <algorithm>
#include <iostream>

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
    psi = psi_b;
}
