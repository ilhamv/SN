#include <cmath>
#include <limits> 
#include <iostream>

#include "Solver.h"
#include "Algorithm.h"


//==============================================================================
// Steady state
//==============================================================================

void source_iteration( const double epsilon,
                       const std::vector<std::shared_ptr<Region>>& mesh,
                       const std::vector<double>& mu, 
                       const std::vector<double>& w,
                       const std::shared_ptr<BC>& BC_left,
                       const std::shared_ptr<BC>& BC_right,
                       AcceleratorDSA& DSA,
                       std::vector<double>& phi,
                       std::vector<double>& rho )
{
    // phi: cell-average scalar flux
    // rho: spectral radius
    const int N = mu.size();
    const int J = phi.size();

    std::vector<double> S(J);        // RHS source (Note: isotropic)
    std::vector<double> psi( N/2 );  // In/out angular flux
    std::vector<double> phi_old(J);
    double psi_avg;                  // Average angular flux
    double error;                    // Maximum relative error
    double rho_num, rho_denom = 1.0;
    int idx;

    do{
        // Set RHS source (Note: isotropic)
        for( int j = 0; j < J; j++ ){
            S[j] = 0.5 * ( mesh[j]->SigmaS() * phi[j] + mesh[j]->Q() );
        }

        // Reset phi
        phi_old = phi; 
        std::fill(phi.begin(), phi.end(), 0.0);

        //======================================================================
        // Forward sweep
        //======================================================================

        // Set BC
        BC_left->set_boundary( psi );
        
        // Sweep
        for( int j = 0; j < J; j++ ){
            idx = 0; // Index for psi, as its size is N/2
            for( int n = 0.5*N; n < N; n++ ){
                psi_avg = psi[idx];
                psi[idx] = ( ( mu[n] - 0.5 * mesh[j]->tau() ) * psi[idx] 
                             + S[j] * mesh[j]->dz() )
                           / ( mu[n] + 0.5 * mesh[j]->tau() );
                psi_avg = (psi_avg + psi[idx]) * w[n];
                phi[j] += psi_avg;
                idx++;
            }
        }
        
        //======================================================================
        // Backward sweep
        //======================================================================

        // Set BC
        BC_right->set_boundary( psi );

        // Sweep
        for( int j = J-1; j >= 0; j-- ){
            for( int n = 0; n < 0.5*N; n++ ){
                psi_avg = psi[n];
                psi[n] = ( ( -mu[n] - 0.5 * mesh[j]->tau() ) * psi[n] 
                           + S[j] * mesh[j]->dz() )
                         / ( -mu[n] + 0.5 * mesh[j]->tau() );
                psi_avg = (psi_avg + psi[n]) * w[n];
                phi[j] += psi_avg;
            }
            phi[j] *= 0.5;
        }
        
        // DSA
        DSA.accelerate( mesh, phi_old, phi );
        
        //======================================================================
        // Maximum relative error and spectral radius estimate
        //======================================================================
       
        error = 0.0;
        for( int j = 0; j < J; j++ ){
            // Now phi_old holds the absolute difference between iterates
            phi_old[j] = std::abs( phi[j] - phi_old[j] );
            double val = phi_old[j] / 
                  ( std::abs(phi[j]) + std::numeric_limits<double>::epsilon() );
            if( val > error ){ error = val; }
        }

        rho_num = norm_2(phi_old);
        rho.push_back( rho_num/rho_denom );
        rho_denom = rho_num;

    } while ( error > ( 1.0 - rho.back() ) * epsilon );
}


//==============================================================================
// Time Depndent
//==============================================================================

void source_iteration_TD( const double epsilon,
                          const std::vector<std::shared_ptr<Region>>& mesh,
                          const std::vector<double>& mu, 
                          const std::vector<double>& w,
                          const std::shared_ptr<BC>& BC_left,
                          const std::shared_ptr<BC>& BC_right,
                          AcceleratorDSA& DSA, 
                          const double speed, const double dt, const int K,
                          const std::vector<std::vector<double>>& psi_initial,
                          std::vector<std::vector<double>>& phi )
{
    // Material total and absorption cross section have been time-augmented.
    // psi_initial: initial cell-edge angular flux [j][n]
    // phi: cell-averaged scalar flux [k][j]
    // K: # of time steps
    const int N = mu.size();
    const int J = phi[0].size();

    // Source iteration tools
    std::vector<double> psi( N/2 );           // In/out angular flux
    std::vector<std::vector<double>> psi_avg; // Cell-average angular flux
    std::vector<std::vector<double>> psi_prv; // Previous cell-average ang flux
    std::vector<double> S(J);                 // Isotropic source
    std::vector<double> phi_old(J);
    double error;                             // Maximum relative error
    double rho;                               // Spectral radius estimate
    double rho_num, rho_denom;
    int idx;
    const double aug = 1.0 / speed / dt;      // Factor in previous time source
    
    // Initial cell-average angular flux
    psi_avg.resize( J, std::vector<double>(N,0.0) );
    psi_prv.resize( J, std::vector<double>(N,0.0) );
    for( int j = 0; j < J; j++ ){
        for( int n = 0; n < N; n++ ){
            psi_avg[j][n] = 0.5 * ( psi_initial[j][n] + psi_initial[j+1][n] );
        }
    }

    // Initial cell-average scalar flux
    for( int j = 0; j < J; j++ ){
        for( int n = 0; n < N; n++ ){
            phi[0][j] += psi_avg[j][n]* w[n];
        }
    }

    // Go over time steps
    for( int k = 1; k < K+1; k++ ){

        //======================================================================
        // Source iteration
        //======================================================================

        rho_denom = 1.0;

        // Set first guess and previous time
        psi_prv = psi_avg;
        phi[k] = phi[k-1];

        int N_iter = 0;

        // Iteration start
        do{

            // Set isotropic source
            for( int j = 0; j < J; j++ ){
                S[j] = 0.5 * ( mesh[j]->SigmaS() * phi[k][j] + mesh[j]->Q() );
            }

            // Reset phi
            phi_old = phi[k];
            std::fill(phi[k].begin(), phi[k].end(), 0.0);

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
                                 + ( S[j] + aug * psi_prv[j][n] )
                                   * mesh[j]->dz() )
                               / ( mu[n] + 0.5 * mesh[j]->tau() );
                    psi_avg[j][n] = 0.5 * (psi_avg[j][n] + psi[idx]);
                    phi[k][j] += psi_avg[j][n] * w[n];
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
                               + ( S[j] + aug * psi_prv[j][n] )
                                   * mesh[j]->dz() )
                             / ( -mu[n] + 0.5 * mesh[j]->tau() );
                    psi_avg[j][n] = 0.5 * (psi_avg[j][n] + psi[n]);
                    phi[k][j] += psi_avg[j][n] * w[n];
                }
            }

            // DSA
            //DSA.accelerate( mesh, phi_old, phi[k] );
            
            //==================================================================
            // Maximum relative error and spectral radius estimate
            //==================================================================
           
            error = 0.0;
            for( int j = 0; j < J; j++ ){
                // Now phi_old holds the absolute difference between iterates
                phi_old[j] = std::abs( phi[k][j] - phi_old[j] );
                double val = phi_old[j] / ( std::abs(phi[k][j]) 
                        + std::numeric_limits<double>::epsilon() );
                if( val > error ){ error = val; }
            }

            rho_num = norm_2(phi_old);
            rho = rho_num/rho_denom;
            rho_denom = rho_num;

            N_iter++;

        } while ( error > ( 1.0 - rho ) * epsilon );
        
        // Some outputs
        std::cout<< "Report for k = " << k << "\n";
        std::cout<< "  Number of iterations: " << N_iter << "\n";
        std::cout<< "  Iteration matrix spectral radius: " << rho << "\n";
    }
}
