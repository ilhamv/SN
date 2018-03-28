#include <cmath>
#include <limits> 

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

void transport_sweep( const double epsilon,
                      const std::vector<std::shared_ptr<Region>>& mesh,
                      const std::vector<double>& mu, 
                      const std::vector<double>& w,
                      const std::shared_ptr<BC>& BC_left,
                      const std::shared_ptr<BC>& BC_right,
                      AcceleratorDSA& DSA, 
                      const double speed, const double dt,
                      std::vector<double>& phi,
                      double& rho,
                      std::vector<std::vector<double>>& psi_prv,
                      std::vector<std::vector<double>>& psi_avg,
                      int& N_iter )
{
    const int N = mu.size();
    const int J = phi.size();

    // Source iteration tools
    std::vector<double> psi( N/2 );           // In/out angular flux
    std::vector<double> S(J);                 // Isotropic source
    std::vector<double> phi_old(J);
    double error;                             // Maximum relative error
    double rho_num, rho_denom = 1.0;
    int idx;

    // Iteration start
    do{

        // Set isotropic source
        for( int j = 0; j < J; j++ ){
            S[j] = 0.5 * ( mesh[j]->SigmaS() * phi[j] + mesh[j]->Q() );
        }

        // Reset phi
        phi_old = phi;
        std::fill(phi.begin(), phi.end(), 0.0);

        //==================================================================
        // Forward sweep
        //==================================================================

        // Set BC
        BC_left->set_boundary( psi );
        
        // Sweep
        for( int j = 0; j < J; j++ ){
            idx = 0; // Index for psi, as its size is N/2
            for( int n = 0.5*N; n < N; n++ ){
                psi_avg[j][n] = psi[idx];
                psi[idx] = ( ( mu[n] - 0.5 * mesh[j]->tau() ) * psi[idx] 
                             + ( S[j] + psi_prv[j][n] )
                               * mesh[j]->dz() )
                           / ( mu[n] + 0.5 * mesh[j]->tau() );
                psi_avg[j][n] = 0.5 * (psi_avg[j][n] + psi[idx]);
                phi[j] += psi_avg[j][n] * w[n];
                idx++;
            }
        }
        
        //==================================================================
        // Backward sweep
        //==================================================================

        // Set BC
        BC_right->set_boundary( psi );

        // Sweep
        for( int j = J-1; j >= 0; j-- ){
            for( int n = 0; n < 0.5*N; n++ ){
                psi_avg[j][n] = psi[n];
                psi[n] = ( ( -mu[n] - 0.5 * mesh[j]->tau() ) * psi[n] 
                           + ( S[j] + psi_prv[j][n] )
                               * mesh[j]->dz() )
                         / ( -mu[n] + 0.5 * mesh[j]->tau() );
                psi_avg[j][n] = 0.5 * (psi_avg[j][n] + psi[n]);
                phi[j] += psi_avg[j][n] * w[n];
            }
        }

        // DSA
        DSA.accelerate( mesh, phi_old, phi );
        
        //==================================================================
        // Maximum relative error and spectral radius estimate
        //==================================================================
       
        error = 0.0;
        for( int j = 0; j < J; j++ ){
            // Now phi_old holds the absolute difference between iterates
            phi_old[j] = std::abs( phi[j] - phi_old[j] );
            double val = phi_old[j] / ( std::abs(phi[j]) 
                    + std::numeric_limits<double>::epsilon() );
            if( val > error ){ error = val; }
        }

        rho_num = norm_2(phi_old);
        rho = rho_num/rho_denom;
        rho_denom = rho_num;

        N_iter++;

    } while ( error > ( 1.0 - rho ) * epsilon );
}
