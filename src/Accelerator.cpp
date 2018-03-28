#include <iostream>
#include "Accelerator.h"


AcceleratorDSA::AcceleratorDSA
        ( const std::vector<std::shared_ptr<Region>>& mesh,
          std::shared_ptr<BC> BC_left,
          std::shared_ptr<BC> BC_right )
{
    //==========================================================================
    // Set the original tridiagonal matrix
    //==========================================================================

    J = mesh.size();
    A.resize(J);
    B.resize(J+1,0.0);
    C.resize(J);

    for( int j = 0; j < J; j++ ){
        A[j] = 0.25 * mesh[j]->SigmaA() * mesh[j]->dz()
               - 1.0 / ( 3.0 * mesh[j]->tau() );
        B[j] = 0.25 * mesh[j]->SigmaA() * mesh[j]->dz()
               + 1.0 / ( 3.0 * mesh[j]->tau() );
    }
    for( int j = 1; j < J+1; j++ ){
        B[j] += 0.25 * mesh[j-1]->SigmaA() * mesh[j-1]->dz()
                + 1.0 / ( 3.0 * mesh[j-1]->tau() );
    }
    C = A;

    B[0] *= 2.0;
    C[0] *= 2.0;
    if( BC_left->type() == "Vacuum" ){ B[0] += 1.0; }
    
    A[J-1] *= 2.0;
    B[J] *= 2.0;
    if( BC_right->type() == "Vacuum" ){ B[J] += 1.0; }

    
    //==========================================================================
    // LU Factorization
    //==========================================================================

    for( int j = 0; j < J; j++ ){
        A[j] /= B[j];
        B[j+1] -= A[j] * C[j];
    }
}


//==============================================================================
// Accelerate
//==============================================================================

void AcceleratorDSA::accelerate(const std::vector<std::shared_ptr<Region>>& mesh,
                                const std::vector<double>& phi_old, 
                                      std::vector<double>& phi )
{
    // RHS
    std::vector<double> f(J+1);
    f[0] = mesh[0]->SigmaS() * mesh[0]->dz() * ( phi[0] - phi_old[0] );
    for( int j = 1; j < J; j++ ){
        f[j] = 0.5 * ( mesh[j-1]->SigmaS() * mesh[j-1]->dz() 
                       * ( phi[j-1] - phi_old[j-1] )
                     + mesh[j]->SigmaS() * mesh[j]->dz() 
                       * ( phi[j] - phi_old[j] ) );
    }
    f[J] = mesh[J-1]->SigmaS() * mesh[J-1]->dz() * ( phi[J-1] - phi_old[J-1] );

    // Forward substitution
    for( int j = 1; j < J+1; j++ ){
        f[j] -= A[j-1] * f[j-1];
    }
    
    // Backward substitution
    f[J] /= B[J];
    for( int j = J-1; j >= 0; j-- ){
        f[j] = ( f[j] - C[j] * f[j+1] ) / B[j];
    }

    // Update phi
    for( int j = 0; j < J; j++ ){
        phi[j] = phi[j] + 0.5 * ( f[j] + f[j+1] );
    }
}
