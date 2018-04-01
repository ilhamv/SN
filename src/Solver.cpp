#include <cmath>
#include <limits> 
#include <iostream>

#include "Solver.h"
#include "Algorithm.h"
#include "Accelerator.h"


//==============================================================================
// Steady state
//==============================================================================

void source_iteration( int& N_iter,
                       const double epsilon,
                       const std::vector<std::shared_ptr<Region>>& mesh,
                       const std::vector<double>& mu, 
                       const std::vector<double>& w,
                       const std::shared_ptr<BC>& BC_left,
                       const std::shared_ptr<BC>& BC_right,
                       std::vector<double>& phi,
                       std::vector<double>& rho,
                       const std::string accelerator_type )
{
    // phi: cell-average scalar flux
    // rho: spectral radius[l]

    // Report mode
    std::cout<< "Mode: Steady-state\n";

    // Tools
    const int N = mu.size();
    const int J = mesh.size();
    std::vector<double> S(J);        // Isotropic source
    std::vector<double> psi( N/2 );  // In/out angular flux
    std::vector<double> phi_old(J);
    double psi_avg;                  // Average angular flux
    double error;                    // Maximum relative error
    double rho_num, rho_denom = 1.0;
    int idx;
    
    // Set accelerator
    std::shared_ptr<Accelerator> accelerator;
    std::cout<< "Accelerator: ";
    if( accelerator_type == "DSA" )
    { 
        std::cout<< "DSA\n";
        accelerator = std::make_shared<AcceleratorDSA>( mesh, BC_left, 
                                                        BC_right );
    } else{
        std::cout<< "OFF\n";
        accelerator = std::make_shared<AcceleratorNONE>();
    }

    // Initialize phi, zero first guess
    phi.resize( J, 0.0 );
   
    //=========================================================================
    // Source iteration
    //=========================================================================

    do{
        // Set RHS isotropic source
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
        
        // Acceleration
        accelerator->accelerate( mesh, phi_old, phi );
        
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
        
    // Some outputs
    N_iter = rho.size();
    std::cout<< "Steady-state source iteration done!\n";
    std::cout<< "Number of iterations: " << N_iter << "\n";
    std::cout<< "Iteration matrix spectral radius: " << rho.back() << "\n";
}


//==============================================================================
// Time Depndent - Implicit
//==============================================================================

void source_iteration_TD_implicit( 
        const double epsilon,
        const std::vector<std::shared_ptr<Region>>& mesh,
        const std::vector<std::shared_ptr<Material>>& material,
        const std::vector<std::shared_ptr<Region>>& region,
        const std::vector<double>& mu, 
        const std::vector<double>& w,
        const std::shared_ptr<BC>& BC_left,
        const std::shared_ptr<BC>& BC_right,
        const double speed, const double dt,
        const std::vector<std::vector<double>>& psi_initial,
        std::vector<std::vector<double>>& phi,
        const std::string accelerator_type,
        const std::vector<double>& time )
{
    // psi_initial: initial cell-edge angular flux [j][n]
    // phi: cell-averaged scalar flux [k][j]

    // Report mode
    std::cout<< "Mode: Time-dependent (Implicit)\n";

    // Source iteration tools
    const int N = mu.size();
    const int J = mesh.size();
    const int K = time.size() - 1;
    std::vector<std::vector<double>> psi_avg; // Cell-average angular flux
    std::vector<std::vector<double>> psi_prv; // Previous cell-average ang flux
    double rho;                               // Spectral radius estimate
    std::vector<double> psi( N/2 );           // In/out angular flux
    std::vector<double> S(J);                 // Isotropic source
    std::vector<double> phi_old(J);
    double error;                             // Maximum relative error
    double rho_num, rho_denom = 1.0;
    int idx;

    // Time augment absorption and total cross sections
    const double aug = 1.0 / speed / dt;
    for( int i = 0; i < material.size(); i++ ){
        material[i]->time_augment( aug );
    }
    for( int i = 0; i < region.size(); i++ ){
        region[i]->reset();
    }

    // Set accelerator
    std::shared_ptr<Accelerator> accelerator;
    std::cout<< "Accelerator: ";
    if( accelerator_type == "DSA" )
    { 
        std::cout<< "DSA\n";
        accelerator = std::make_shared<AcceleratorDSA>( mesh, BC_left, 
                                                        BC_right );
    } else{
        std::cout<< "OFF\n";
        accelerator = std::make_shared<AcceleratorNONE>();
    }


    //=========================================================================
    // Initial conditions
    //=========================================================================

    // Initialize
    psi_avg.resize( J, std::vector<double>(N,0.0) );
    psi_prv.resize( J, std::vector<double>(N,0.0) );
    phi.resize( K+1, std::vector<double>(J,0.0) );
    
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


    //=========================================================================
    // March in time
    //=========================================================================
    
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


        //======================================================================
        // Source Iteration - Implicit time dependent
        //======================================================================
        
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
                                 + ( S[j] + psi_prv[j][n] )
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
                               + ( S[j] + psi_prv[j][n] )
                                   * mesh[j]->dz() )
                             / ( -mu[n] + 0.5 * mesh[j]->tau() );
                    psi_avg[j][n] = 0.5 * (psi_avg[j][n] + psi[n]);
                    phi[k][j] += psi_avg[j][n] * w[n];
                }
            }

            // Accelerate
            accelerator->accelerate( mesh, phi_old, phi[k] );
            
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
        std::cout<< "Report for k = " << k << " ("<< time[k] << " s)\n";
        std::cout<< "  Number of iterations: " << N_iter << "\n";
        std::cout<< "  Iteration matrix spectral radius: " << rho << "\n";
    }
        
    // Revert time augment
    for( int i = 0; i < material.size(); i++ ){
        material[i]->revert_augment( aug );
    }
    for( int i = 0; i < region.size(); i++ ){
        region[i]->reset();
    }
}

/*
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
        region[i]->reset();
    }
    
    AcceleratorDSA DSA_aug( mesh, BC_left, BC_right );
    
    // Revert time augment
    for( int i = 0; i < material.size(); i++ ){
        material[i]->revert_augment( aug );
    }
    for( int i = 0; i < region.size(); i++ ){
        region[i]->reset();
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

        // Set previous time
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

            // Set up first guess
            phi[k] = phi[k-1];

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
                region[i]->reset();
            }
            
            // Source iteration
            transport_sweep( epsilon, mesh, mu, w, BC_left, BC_right, DSA_aug, 
                             speed, dt, phi[k], rho_add, S_anis,psi_avg,N_iter);

            // Revert time augment
            for( int i = 0; i < material.size(); i++ ){
                material[i]->revert_augment( aug );
            }
            for( int i = 0; i < region.size(); i++ ){
                region[i]->reset();
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
*/


//==============================================================================
// Time Depndent - Multiple Balance
//==============================================================================

void source_iteration_TD_MB( 
        const double epsilon,
        const std::vector<std::shared_ptr<Region>>& mesh,
        const std::vector<std::shared_ptr<Material>>& material,
        const std::vector<std::shared_ptr<Region>>& region,
        const std::vector<double>& mu, 
        const std::vector<double>& w,
        const std::shared_ptr<BC>& BC_left,
        const std::shared_ptr<BC>& BC_right,
        const double speed, const double dt,
        const std::vector<std::vector<double>>& psi_initial,
        std::vector<std::vector<double>>& phi,
        const std::string accelerator_type,
        const std::vector<double>& time )
{;}/*
    // psi_initial: initial cell-edge angular flux [j][n]
    // phi: cell-averaged scalar flux [k][j]

    // Report mode
    std::cout<< "Mode: Time-dependent (Multiple-Balance)\n";

    // Source iteration tools
    const int N = mu.size();
    const int J = mesh.size();
    const int K = time.size() - 1;
    std::vector<std::vector<double>> psi_avg; // Cell-average angular flux
    std::vector<std::vector<double>> psi_prv; // Previous cell-average ang flux
    std::vector<std::vector<double>> psi_prv; // Previous cell-average ang flux
    double rho;                               // Spectral radius estimate
    std::vector<double> psi( N/2 );           // In/out angular flux
    std::vector<double> phi_old(J);
    std::vector<double> J(J+1);               // Cell-edge current
    std::vector<double> J_old(J+1);
    std::vector<std::vector<double>> S;       // RHS source
    double error;                             // Maximum relative error
    double rho_num, rho_denom = 1.0;
    int idx;
    double Sb, A, B, C;

    // Set accelerator
    std::shared_ptr<Accelerator> accelerator;
    std::cout<< "Accelerator: ";
    if( accelerator_type == "DSA" )
    { 
        std::cout<< "DSA for MB is not available yet --> OFF\n";
        accelerator = std::make_shared<AcceleratorNONE>();
    } else{
        std::cout<< "OFF\n";
        accelerator = std::make_shared<AcceleratorNONE>();
    }


    //=========================================================================
    // Initial conditions
    //=========================================================================

    // Initialize
    psi_avg.resize( J, std::vector<double>(N,0.0) );
    psi_prv.resize( J, std::vector<double>(N,0.0) );
    phi.resize( K+1, std::vector<double>(J,0.0) );
    S.resize( J-1, std::vector<double>(N,0.0) );
    
    // Initial cell-average
    for( int j = 0; j < J; j++ ){
        for( int n = 0; n < N; n++ ){
            // Angular flux
            psi_avg[j][n] = 0.5 * ( psi_initial[j][n] + psi_initial[j+1][n] );
            
            // Scalar flux
            phi[0][j] += psi_avg[j][n]* w[n];
        }
    }

    // Initial cell-edge current
    for( int j = 0; j < J+1; j++ ){
        for( int n = 0; n < N; n++ ){
            J[j] += mu[n] * psi_initial[j][n]* w[n];
        }
    }


    //=========================================================================
    // March in time
    //=========================================================================
    
    for( int k = 1; k < K+1; k++ ){
        // Iteration counter
        int N_iter = 0;

        // Set first guess and previous time
        psi_prv = psi_avg;
        phi[k] = phi[k-1];
        J_old = J;

        
        //======================================================================
        // Source Iteration - Implicit time dependent
        //======================================================================
        
        do{
            // Set RHS source
            for( int j = 0; j < J-1; j++ ){
                for( int n = 0; n < N; n++ ){
                    S[j][n] = ( mesh[j+1]->dz * ( 1.0 / ( speed * dt ) + 0.5 * ( mesh[j+1]->SigmaT() + mesh[j+1]->SigmaA() ) ) + mu[n] ) *  mesh[j+1]->SigmaS() * phi[k][j+1]
                              + ( mesh[j]->dz * ( 1.0 / ( speed * dt ) + 0.5 * ( mesh[j]->SigmaT() + mesh[j]->SigmaA() ) ) - mu[n] ) *  mesh[j]->SigmaS() * phi[k][j]
                              + 0.5 * ( mesh[j+1]->SigmaS() * J[j+2] + ( mesh[j]->SigmaS() - mesh[j+1]->SigmaS() ) J[j+1] - mesh[j]->SigmaS() * J[j] )
                              + ( mesh[j+1]->dz * ( 1.0 / ( speed * dt ) + 0.5 * mesh[j+1]->SigmaA() ) + mu[n] ) * mesh[j+1]->Q()
                              + ( mesh[j]->dz * ( 1.0 / ( speed * dt ) + 0.5 * mesh[j]->SigmaA() ) - mu[n] ) * mesh[j]->Q()
                              + 2.0 / ( speed * speed * dt * dt ) * ( mesh[j+1]->dz() * psi_prv[j+1][n] + mesh[j]->dz() * psi_prv[j][n] );

            }

            // Reset phi and J
            phi_old = phi[k];
            J_old = J;

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
                }
            }

            // Accelerate
            accelerator->accelerate( mesh, phi_old, phi[k] );

            //==================================================================
            // Update phi and J
            //==================================================================
            
            for( int j = 0; j < J-1; j++ ){
                phi[k][j] = 0.0;
                J[j] = 0.0;
                for( int n = 0; n < N; n++ ){
                    phi[k][j] += psi_avg[j][n] * w[n];
                    J[j] = mu[n] * psi_avg[j][n] * w[n];
                }
            }

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
        std::cout<< "Report for k = " << k << " ("<< time[k] << " s)\n";
        std::cout<< "  Number of iterations: " << N_iter << "\n";
        std::cout<< "  Iteration matrix spectral radius: " << rho << "\n";
    }
    }
}
*/
