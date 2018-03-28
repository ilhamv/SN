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
    
    
    // Set up DSA
    AcceleratorDSA DSA( mesh, BC_left, BC_right );

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
                        const std::vector<std::shared_ptr<Material>>& material,
                          const std::vector<std::shared_ptr<Region>>& region,
                          const std::vector<double>& mu, 
                          const std::vector<double>& w,
                          const std::shared_ptr<BC>& BC_left,
                          const std::shared_ptr<BC>& BC_right,
                          const double speed, const double dt, const int K,
                          const std::vector<std::vector<double>>& psi_initial,
                          std::vector<std::vector<double>>& phi )
{
    // psi_initial: initial cell-edge angular flux [j][n]
    // phi: cell-averaged scalar flux [k][j]
    // K: # of time steps
    const int N = mu.size();
    const int J = phi[0].size();

    // Source iteration tools
    std::vector<std::vector<double>> psi_avg; // Cell-average angular flux
    std::vector<std::vector<double>> psi_prv; // Previous cell-average ang flux
    double rho;                               // Spectral radius estimate
    const double aug = 1.0 / speed / dt;      // Factor in previous time source
    
    // Time augment absorption and total cross sections
    for( int i = 0; i < material.size(); i++ ){
        material[i]->time_augment( aug );
    }
    for( int i = 0; i < region.size(); i++ ){
        region[i]->reset_tau();
    }

    // Set Accelerator
    AcceleratorDSA DSA( mesh, BC_left, BC_right );
    
    //=========================================================================
    // Set initial conditions
    //=========================================================================

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
        // Iteration counter
        int N_iter = 0;

        // Set first guess and previous time
        psi_prv = psi_avg;
        phi[k] = phi[k-1];

        // Set up psi_prv as the non isotropic source
        for( int j = 0; j < J; j++ ){
            for( int n = 0; n < N; n++ ){
                psi_prv[j][n] *= aug;
            }
        }

        // Source iteration
        transport_sweep( epsilon, mesh, mu, w, BC_left, BC_right, DSA, 
                         speed, dt, phi[k], rho, psi_prv, psi_avg, N_iter );
        
        // Some outputs
        std::cout<< "Report for k = " << k << "\n";
        std::cout<< "  Number of iterations: " << N_iter << "\n";
        std::cout<< "  Iteration matrix spectral radius: " << rho << "\n";
    }

    // Revert time augment
    for( int i = 0; i < material.size(); i++ ){
        material[i]->revert_augment( aug );
    }
    for( int i = 0; i < region.size(); i++ ){
        region[i]->reset_tau();
    }
}


//==============================================================================
// Time Depndent - Multiple Balance
//==============================================================================

void source_iteration_MB( const double epsilon,
                          const std::vector<std::shared_ptr<Region>>& mesh,
                        const std::vector<std::shared_ptr<Material>>& material,
                          const std::vector<std::shared_ptr<Region>>& region,
                          const std::vector<double>& mu, 
                          const std::vector<double>& w,
                          const std::shared_ptr<BC>& BC_left,
                          const std::shared_ptr<BC>& BC_right,
                          const double speed, const double dt, const int K,
                          const std::vector<std::vector<double>>& psi_initial,
                          std::vector<std::vector<double>>& phi )
{
    // psi_initial: initial cell-edge angular flux [j][n]
    // phi: cell-averaged scalar flux [k][j]
    // K: # of time steps
    const int N = mu.size();
    const int J = phi[0].size();

    // Source iteration tools
    std::vector<std::vector<double>> psi_avg; // Cell-average angular flux
    std::vector<std::vector<double>> psi_prv; // Previous cell-average ang flux
    std::vector<std::vector<double>> psi_nxt; // Next cell-average angular flux
    std::vector<std::vector<double>> S_anis;  // Anisotropic source
    double rho_ori;                           // Spectral radius
    double rho_add;                           // Spectral radius
    const double aug = 2.0 / speed / dt;      // Factor in previous time source
    const double aug_half = 1.0 / speed / dt; // Factor in previous time source
    double error;                             // Maximum relative error
   
    //=========================================================================
    // Set accelerator for each transport
    //=========================================================================

    AcceleratorDSA DSA( mesh, BC_left, BC_right );

    // Time augment absorption and total cross sections
    for( int i = 0; i < material.size(); i++ ){
        material[i]->time_augment( aug );
    }
    for( int i = 0; i < region.size(); i++ ){
        region[i]->reset_tau();
    }
    
    AcceleratorDSA DSA_aug( mesh, BC_left, BC_right );
    
    // Revert time augment
    for( int i = 0; i < material.size(); i++ ){
        material[i]->revert_augment( aug );
    }
    for( int i = 0; i < region.size(); i++ ){
        region[i]->reset_tau();
    }
    
    //=========================================================================
    // Set initial conditions
    //=========================================================================

    // Initialize
    psi_avg.resize( J, std::vector<double>(N,0.0) );
    psi_prv.resize( J, std::vector<double>(N,0.0) );
    psi_nxt.resize( J, std::vector<double>(N,0.0) );
    S_anis.resize( J, std::vector<double>(N,0.0) );

    // Initial cell-average angular flux
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
        // Iteration counter
        int N_iter = 0;
        int N_iter_out = 0;

        // Set first guesses and previous time
        phi[k] = phi[k-1];
        psi_prv = psi_avg;

        do{
            // Reset next cell-edge angular flux
            psi_nxt = psi_avg;

            //==================================================================
            // Original transport
            //==================================================================

            // Time statuses:
            //   psi_avg: --> time-average (for additional transport)
            //   psi_prv: previous (remain unchanged)
            //   psi_nxt: next (previous estimate)
            //   phi[k]:  --> time-average (not used)

            // Set up anisotropic source
            for( int j = 0; j < J; j++ ){
                for( int n = 0; n < N; n++ ){
                    S_anis[j][n] = aug_half * ( psi_prv[j][n] - psi_nxt[j][n] );
                }
            }

            // Source iteration
            transport_sweep( epsilon, mesh, mu, w, BC_left, BC_right, DSA, 
                             speed, dt, phi[k], rho_ori, S_anis,psi_avg,N_iter);


            //==================================================================
            // Additional transport
            //==================================================================
            
            // Time statuses:
            //   psi_avg: time-average (RHS source) --> next (new estimate)
            //   psi_nxt: next (previous estimate)
            //   phi[k]:  --> next (result we want)

            // Set up anisotropic source
            for( int j = 0; j < J; j++ ){
                for( int n = 0; n < N; n++ ){
                    S_anis[j][n] = aug * psi_avg[j][n];
                }
            }

            // Time augment absorption and total cross sections
            for( int i = 0; i < material.size(); i++ ){
                material[i]->time_augment( aug );
            }
            for( int i = 0; i < region.size(); i++ ){
                region[i]->reset_tau();
            }
            
            // Source iteration
            transport_sweep( epsilon, mesh, mu, w, BC_left, BC_right, DSA_aug, 
                             speed, dt, phi[k], rho_add, S_anis,psi_avg,N_iter);
            
            // Revert time augment
            for( int i = 0; i < material.size(); i++ ){
                material[i]->revert_augment( aug );
            }
            for( int i = 0; i < region.size(); i++ ){
                region[i]->reset_tau();
            }
            
            
            //==================================================================
            // Maximum relative error
            //==================================================================
           
            error = 0.0;
            for( int j = 0; j < J; j++ ){
                for( int n = 0; n < N; n++ ){
                    // Now psi_nxt holds the absolute difference between iterate
                    psi_nxt[j][n] = std::abs( psi_avg[j][n] - psi_nxt[j][n] );
                    double val = psi_nxt[j][n] / ( std::abs( psi_avg[j][n] ) 
                                 + std::numeric_limits<double>::epsilon() );
                    if( val > error ){ error = val; }
                }
            }

            N_iter_out++;

        } while ( error > epsilon );

        // Some outputs
        std::cout<< "Report for k = " << k << "\n";
        std::cout<< "  Number of iterations: " << N_iter << "\n";
        std::cout<< "  Number of outer iterations: " << N_iter_out << "\n";
        std::cout<< "  Original transport spectral radius: " << rho_ori << "\n";
        std::cout<< "  Additional transport spectral radius: " << rho_add<<"\n";
    }

}
