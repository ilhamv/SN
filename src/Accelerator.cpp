#include <iostream>
#include "Accelerator.h"


//==============================================================================
// Set Consistent DSA
//==============================================================================

AcceleratorDSA::AcceleratorDSA
        ( const std::vector<std::shared_ptr<Region>>& mesh,
          std::shared_ptr<BC> BC_left,
          std::shared_ptr<BC> BC_right ): Accelerator("DSA")
{
    //==========================================================================
    // Set the original tridiagonal matrix
    //==========================================================================

    J = mesh.size();
    A.resize(J);
    B.resize(J+1);
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
    if( BC_left->type() != "Reflective" ){ B[0] += 1.0; }
    
    A[J-1] *= 2.0;
    B[J] *= 2.0;
    if( BC_right->type() != "Reflective" ){ B[J] += 1.0; }

    
    //==========================================================================
    // LU Factorization
    //==========================================================================

    for( int j = 0; j < J; j++ ){
        A[j] /= B[j];
        B[j+1] -= A[j] * C[j];
    }
}


//==============================================================================
// Set Inconsistent DSA
//==============================================================================

AcceleratorIDSA::AcceleratorIDSA
        ( const std::vector<std::shared_ptr<Region>>& mesh,
          std::shared_ptr<BC> BC_left,
          std::shared_ptr<BC> BC_right ): Accelerator("IDSA")
{
    //==========================================================================
    // Set the original tridiagonal matrix
    //==========================================================================

    J = mesh.size();
    A.resize(J-1);
    B.resize(J);
    C.resize(J-1);

    double SigmaT, h, val;

    for( int j = 0; j < J-1; j++ ){
        SigmaT = ( mesh[j]->tau() + mesh[j+1]->tau() )
                 / ( mesh[j]->dz() + mesh[j+1]->dz() );
        h = 0.5 * ( mesh[j]->dz() + mesh[j+1]->dz() );
        val = 1.0 / ( 3.0 * SigmaT * h );
        A[j] = -val;
        B[j] = val;
    }
    C = A;
    for( int j = 1; j < J; j++ ){
        SigmaT = ( mesh[j]->tau() + mesh[j-1]->tau() )
                 / ( mesh[j]->dz() + mesh[j-1]->dz() );
        h = 0.5 * ( mesh[j]->dz() + mesh[j-1]->dz() );
        val = 1.0 / ( 3.0 * SigmaT * h );
        B[j] += val;
    }
    for( int j = 0; j < J; j++ ){
        B[j] += mesh[j]->SigmaA() * mesh[j]->dz();
    }

    // Left BC
    SigmaT = mesh[0]->SigmaT();
    h = 0.5 * mesh[0]->dz();
    val = 1.0 / ( 3.0 * SigmaT * h );
    if( BC_left->type() != "Reflective" ){ 
        B[0] += val * ( 1.0 - 4.0 / ( 4.0 + 3.0 * mesh[0]->tau() ) );
    }
   
    // Right BC
    SigmaT = mesh[J-1]->SigmaT();
    h = 0.5 * mesh[J-1]->dz();
    val = 1.0 / ( 3.0 * SigmaT * h );
    if( BC_right->type() != "Reflective" ){ 
        B[J-1] += val * ( 1.0 - 4.0 / ( 4.0 + 3.0 * mesh[J-1]->tau() ) );
    }

    
    //==========================================================================
    // LU Factorization
    //==========================================================================

    for( int j = 0; j < J-1; j++ ){
        A[j] /= B[j];
        B[j+1] -= A[j] * C[j];
    }
}


//==============================================================================
// No acceleration
//==============================================================================

void AcceleratorNONE::accelerate(const std::vector<std::shared_ptr<Region>>& mesh,
                                 const std::vector<double>& phi_old, 
                                       std::vector<double>& phi ) { return; }

//==============================================================================
// Consistent DSA acceleration
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

//==============================================================================
// Inconsistent DSA acceleration
//==============================================================================

void AcceleratorIDSA::accelerate(const std::vector<std::shared_ptr<Region>>& mesh,
                                const std::vector<double>& phi_old, 
                                      std::vector<double>& phi )
{
    // RHS
    std::vector<double> f(J);
    for( int j = 0; j < J; j++ ){
        f[j] = mesh[j]->SigmaS() * ( phi[j] - phi_old[j] ) * mesh[j]->dz();
    }

    // Forward substitution
    for( int j = 1; j < J; j++ ){
        f[j] -= A[j-1] * f[j-1];
    }
    
    // Backward substitution
    f[J-1] /= B[J-1];
    for( int j = J-2; j >= 0; j-- ){
        f[j] = ( f[j] - C[j] * f[j+1] ) / B[j];
    }

    // Update phi
    for( int j = 0; j < J; j++ ){
        phi[j] = phi[j] + f[j];
    }
}
