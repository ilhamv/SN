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
                       const std::vector<std::shared_ptr<Region>>& region,
                       const std::vector<double>& mu, 
                       const std::vector<double>& w,
                       const std::shared_ptr<BC>& BC_left,
                       const std::shared_ptr<BC>& BC_right,
                       std::vector<double>& phi,
                       std::vector<std::vector<double>>& psi,
                       std::vector<double>& rho,
                       const std::string space_method,
                       const std::string accelerator_type )
{
    // phi: cell-average scalar flux
    // psi: cell-edge angular flux
    // rho: spectral radius

    // Report mode
    std::cout<< "Mode: Steady-state\n";

    // Space method, and set region's alpha
    std::cout<<"Space discretization: "<<space_method<<"\n";
    for( int i = 0; i < region.size(); i++ ){
        region[i]->set_alpha( mu, space_method );
    }

    // Tools
    const int N = mu.size();
    const int J = mesh.size();
    double              S;             // Isotropic source
    std::vector<double> psi_io( N/2 ); // In/out angular flux
    std::vector<double> phi_old(J);
    double              error;         // Maximum relative error
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
    } else if ( accelerator_type == "IDSA" ){
        std::cout<< "IDSA\n";
        accelerator = std::make_shared<AcceleratorIDSA>( mesh, BC_left, 
                                                         BC_right );
    } else{
        std::cout<< "OFF\n";
        accelerator = std::make_shared<AcceleratorNONE>();
    }

    // Initialize phi and psi (zero first guess)
    phi.resize( J, 0.0 );
    psi.resize( J+1, std::vector<double>(N) );
   
    //=========================================================================
    // Source iteration
    //=========================================================================

    do{
        //======================================================================
        // Forward sweep
        //======================================================================

        // Set BC
        BC_left->set_boundary( psi_io );
        
        // Space sweep
        for( int j = 0; j < J; j++ ){
            idx = 0; // Index for psi_io, as its size is N/2

            // Set isotropic source
            S = 0.5 * ( mesh[j]->SigmaS() * phi[j] + mesh[j]->Q() );

            // Direction sweep
            for( int n = 0.5*N; n < N; n++ ){
                psi[j][n]   = psi_io[idx];
                psi_io[idx] = ( ( mu[n] - ( 1.0 - mesh[j]->alpha(n) ) * 0.5 
                                          * mesh[j]->tau() ) * psi_io[idx] 
                                + S * mesh[j]->dz() )
                              / ( mu[n] + ( 1.0 + mesh[j]->alpha(n) ) * 0.5 
                                          * mesh[j]->tau() );
                psi[j+1][n] = psi_io[idx]; idx++;
            }
        }
        
        //======================================================================
        // Backward sweep
        //======================================================================

        // Set BC
        BC_right->set_boundary( psi_io );

        // Space sweep
        for( int j = J-1; j >= 0; j-- ){
            // Set isotropic source
            S = 0.5 * ( mesh[j]->SigmaS() * phi[j] + mesh[j]->Q() );

            // Direction sweep
            for( int n = 0; n < 0.5*N; n++ ){
                psi[j+1][n] = psi_io[n];
                psi_io[n]   = ( ( -mu[n] - ( 1.0 + mesh[j]->alpha(n) ) * 0.5 
                                           * mesh[j]->tau() ) * psi_io[n] 
                                + S * mesh[j]->dz() )
                              / ( -mu[n] + ( 1.0 - mesh[j]->alpha(n) ) * 0.5 
                                         * mesh[j]->tau() );
                psi[j][n]   = psi_io[n];
            }
        }
        
        //======================================================================
        // Update flux
        //======================================================================
        
        // Store old phi
        phi_old = phi; 

        // Update phi
        for( int j = 0; j < J; j++ ){
            phi[j] = 0.0;
            for( int n = 0; n < N; n++ ){
                phi[j] += ( ( 1.0 - mesh[j]->alpha(n) ) * psi[j][n] 
                            +( 1.0 + mesh[j]->alpha(n) ) * psi[j+1][n] ) * w[n];
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

    } while ( error > ( rho.back() / ( 1.0 - rho.back() ) ) * epsilon );
        
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

        } while ( error > epsilon );
        
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
{
    // psi_initial: initial cell-edge angular flux [j][n]
    // phi: cell-averaged scalar flux [k][j]

    // Report mode
    std::cout<< "Mode: Time-dependent (Multiple-Balance)\n";

    // Source iteration tools
    const int N = mu.size();
    const int J = mesh.size();
    const int K = time.size() - 1;
    std::vector<std::vector<double>> psi_avg; // Cell-average angular flux
    std::vector<std::vector<double>> psi_edg; // Cell-edge angular flux
    std::vector<std::vector<double>> psi_prv; // Previous cell-average ang flux
    double rho;                               // Spectral radius estimate
    std::vector<double> psi_b( N/2 );         // Boundary angular flux
    double              psi_b_k;              // Boundary angular flux (t-avg)
    std::vector<double> phi_old(J);
    std::vector<double> phi_k(J);
    std::vector<double> Jc(J+1,0);            // Cell-edge current
    std::vector<std::vector<double>> S;       // RHS source
    std::vector<double>              Sb(N);   // BC source
    double error;                             // Maximum relative error
    double rho_num, rho_denom = 1.0;
    int idx;

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
    // Some constants
    //=========================================================================

    const double gamma    = 1.0 / speed/ dt;
    const double gamma_sq = gamma * gamma;
    std::vector<std::vector<double>> eta, xi, alpha, beta, A, B, C;
    std::vector<double> Ab(N), Bb(N);

    eta.resize( J, std::vector<double>(N,0.0) );
    xi.resize( J, std::vector<double>(N,0.0) );
    alpha.resize( J, std::vector<double>(N,0.0) );
    beta.resize( J, std::vector<double>(N,0.0) );
    A.resize( J, std::vector<double>(N,0.0) );
    B.resize( J, std::vector<double>(N,0.0) );
    C.resize( J, std::vector<double>(N,0.0) );

    for( int j = 0; j < J; j++ ){
        for( int n = 0; n < N; n++ ){
            eta[j][n] = mesh[j]->tau() + 2.0 * mu[n];
            xi[j][n]  = mesh[j]->tau() - 2.0 * mu[n];
            alpha[j][n] = 0.5 * mesh[j]->SigmaT() + gamma 
                          + mu[n] / mesh[j]->dz();
            beta[j][n]  = 0.5 * mesh[j]->SigmaT() + gamma 
                          - mu[n] / mesh[j]->dz();
        }
    }
    
    for( int j = 0; j < J-1; j++ ){
        for( int n = 0; n < N; n++ ){
            A[j][n] = eta[j+1][n] * alpha[j+1][n] + mesh[j+1]->dz() * gamma_sq;
            B[j][n] = eta[j+1][n] * beta[j+1][n]  + mesh[j+1]->dz() * gamma_sq
                    + xi[j][n]    * alpha[j][n]   + mesh[j]->dz()   * gamma_sq;
            C[j][n] = xi[j][n]    * beta[j][n]    + mesh[j]->dz()   * gamma_sq;
        }
    }

    for( int n = 0.5*N; n < N; n++ ){
        Ab[n] = eta[0][n] * alpha[0][n] + mesh[0]->dz() * gamma_sq;
        Bb[n] = eta[0][n] * beta[0][n]  + mesh[0]->dz() * gamma_sq;
    }

    for( int n = 0; n < 0.5*N; n++ ){
        Ab[n] = xi[J-1][n] * alpha[J-1][n] + mesh[J-1]->dz() * gamma_sq;
        Bb[n] = xi[J-1][n] * beta[J-1][n] + mesh[J-1]->dz() * gamma_sq;
    }

    //=========================================================================
    // Initial conditions
    //=========================================================================

    // Initialize
    psi_avg.resize( J, std::vector<double>(N,0.0) );
    psi_prv.resize( J, std::vector<double>(N,0.0) );
    psi_edg.resize( J+1, std::vector<double>(N,0.0) );
    phi.resize( K+1, std::vector<double>(J,0.0) );
    S.resize( J-1, std::vector<double>(N,0.0) );
    
    // Initial cell-average
    for( int j = 0; j < J; j++ ){
        phi[0][j] = 0.0;
        for( int n = 0; n < N; n++ ){
            // Angular flux
            psi_avg[j][n] = 0.5 * ( psi_initial[j][n] + psi_initial[j+1][n] );
            
            // Scalar flux
            phi[0][j] += psi_avg[j][n]* w[n];
        }
    }

    // Initial cell-edge current
    psi_edg = psi_initial;
    for( int j = 0; j < J+1; j++ ){
        Jc[j] = 0.0;
        for( int n = 0; n < N; n++ ){
            Jc[j] += mu[n] * psi_edg[j][n]* w[n];
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

        //======================================================================
        // Source Iteration - Implicit time dependent
        //======================================================================
        
        do{
            // Set cell-time-average scalar flux --> phi_old
            for( int j = 0; j < J; j++ ){
                phi_k[j] = 0.5 / gamma * 
                             ( 1.0 / mesh[j]->dz() * ( Jc[j+1] - Jc[j] ) 
                               + ( mesh[j]->SigmaA() + 2.0 * gamma ) * phi[k][j]
                               - mesh[j]->Q() );
            }

            // Set RHS source
            for( int j = 0; j < J-1; j++ ){
                for( int n = 0; n < N; n++ ){
                    S[j][n] = gamma * ( mesh[j+1]->dz() * ( mesh[j+1]->SigmaS() 
                                                            * phi_k[j+1] 
                                                            + mesh[j+1]->Q() )
                                        + mesh[j]->dz() * ( mesh[j]->SigmaS() 
                                                            * phi_k[j] 
                                                            + mesh[j]->Q() ) )
                              + 0.5 * ( eta[j+1][n] * ( mesh[j+1]->SigmaS() 
                                                        * phi[k][j+1] 
                                                        + mesh[j+1]->Q() )
                                        + xi[j][n] * ( mesh[j]->SigmaS() 
                                                        * phi[k][j] 
                                                        + mesh[j]->Q() ) );
                }
            }

            // Set BC source
            for( int n = 0.5*N; n < N; n++ ){
                Sb[n] = gamma * mesh[0]->dz() * ( mesh[0]->SigmaS() 
                                                  * phi_k[0] + mesh[0]->Q() )
                          + 0.5 * eta[0][n] * ( mesh[0]->SigmaS() 
                                                * phi[k][0] + mesh[0]->Q() ) ;
            }

            for( int n = 0; n < 0.5*N; n++ ){
                Sb[n] = gamma * mesh[J-1]->dz() * ( mesh[J-1]->SigmaS() 
                                                    * phi_k[J-1] 
                                                    + mesh[J-1]->Q() )
                          + 0.5 * xi[J-1][n] * ( mesh[J-1]->SigmaS() 
                                                 * phi[k][J-1] 
                                                 + mesh[J-1]->Q() ) ;
            }

            // Set old phi for error calculation
            phi_old = phi[k];

            //==================================================================
            // Begin forward sweep
            //==================================================================
            
            // Set BC
            BC_left->set_boundary( psi_b );
            
            idx = 0;
            for( int n = 0.5*N; n < N; n++ ){
                // Store BC
                psi_edg[0][n] = psi_b[idx];

                // Set time-average BC
                if( BC_left->type() == "Reflective" ){
                    psi_b_k = 1.0 / ( 4.0 * mu[n] * gamma )
                        * (-(( xi[0][n] * beta[0][n] + mesh[0]->dz() * gamma_sq)
                             * psi_edg[1][N-1-n] 
                            +( xi[0][n] * alpha[0][n]+ mesh[0]->dz() * gamma_sq)
                             * psi_edg[0][N-1-n] )
                           + gamma * mesh[0]->dz() * ( mesh[0]->SigmaS() 
                                                       * phi_k[0]
                                                       + mesh[0]->Q() )
                           + 0.5 * xi[0][n] * ( mesh[0]->SigmaS() * phi_old[0]
                                                + mesh[0]->Q() )
                           + 2.0 * mesh[0]->dz() * gamma_sq * psi_prv[0][N-1-n]
                           );
                } else{ psi_b_k = psi_b[idx]; }

                psi_edg[1][n] = ( Sb[n] - Bb[n] * psi_edg[0][n] 
                                  + 2.0 * mesh[0]->dz() * gamma_sq 
                                    * psi_prv[0][n]
                                  + 4.0 * mu[n] * gamma * psi_b_k ) / Ab[n];

                idx++;
            }

            //==================================================================
            // Continue forward sweep
            //==================================================================

            for( int j = 0; j < J-1; j++ ){
                for( int n = 0.5*N; n < N; n++ ){
                    psi_edg[j+2][n] = ( S[j][n] 
                                        - B[j][n] * psi_edg[j+1][n] 
                                        - C[j][n] * psi_edg[j][n] 
                                        + 2.0 * gamma_sq
                                          * ( mesh[j+1]->dz() * psi_prv[j+1][n]
                                              + mesh[j]->dz() * psi_prv[j][n] ))
                                      / A[j][n];
                }
            }
            
            //==================================================================
            // Begin backward sweep
            //==================================================================
            
            // Set BC
            BC_right->set_boundary( psi_b );
            
            for( int n = 0; n < 0.5*N; n++ ){
                // Store BC
                psi_edg[J][n] = psi_b[n];
                
                // Set time-average BC
                if( BC_right->type() == "Reflective" ){
                    psi_b_k = - 1.0 / ( 4.0 * mu[n] * gamma )
                        * ( -( eta[J-1][n] * beta[J-1][n] + mesh[J-1]->dz() 
                                                            * gamma_sq ) 
                             * psi_edg[J][N-1-n]
                            -( eta[J-1][n] * alpha[J-1][n] + mesh[J-1]->dz() 
                                                             * gamma_sq ) 
                             * psi_edg[J-1][N-1-n]
                           + gamma * mesh[J-1]->dz() * ( mesh[J-1]->SigmaS() 
                                                       * phi_k[J-1]
                                                       + mesh[J-1]->Q() )
                           + 0.5 * eta[J-1][n] * ( mesh[J-1]->SigmaS() 
                                                   * phi_old[J-1] 
                                                   + mesh[J-1]->Q() )
                           + 2.0 * mesh[J-1]->dz() * gamma_sq 
                             * psi_prv[J-1][N-1-n] );
                } else{ psi_b_k = psi_b[n]; }

                psi_edg[J-1][n] = ( Sb[n] - Ab[n] * psi_edg[J][n]
                                    + 2.0 * mesh[J-1]->dz() * gamma_sq 
                                      * psi_prv[J-1][n]
                                    - 4.0 * mu[n] * gamma * psi_b_k ) / Bb[n];
            }

            //==================================================================
            // Continue backward sweep
            //==================================================================

            for( int j = J-2; j >= 0; j-- ){
                for( int n = 0; n < 0.5*N; n++ ){
                    psi_edg[j][n] = ( S[j][n] 
                                      - A[j][n] * psi_edg[j+2][n]
                                      - B[j][n] * psi_edg[j+1][n] 
                                      + 2.0 * gamma_sq
                                        * ( mesh[j+1]->dz() * psi_prv[j+1][n]
                                            + mesh[j]->dz() * psi_prv[j][n] ))
                                      / C[j][n];
                }
            }

            //==================================================================
            // Update cell-edge current, cell-average angular and scalar flux
            //==================================================================
            
            for( int j = 0; j < J+1; j++ ){
                Jc[j] = 0;
                for( int n = 0; n < N; n++ ){
                    Jc[j] += mu[n] * psi_edg[j][n]* w[n];
                }
            }
            
            for( int j = 0; j < J; j++ ){
                phi[k][j] = 0.0;
                for( int n = 0; n < N; n++ ){
                    psi_avg[j][n] = 0.5 * ( psi_edg[j][n] + psi_edg[j+1][n] );
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

        } while ( error > epsilon );
        
        // Some outputs
        std::cout<< "Report for k = " << k << " ("<< time[k] << " s)\n";
        std::cout<< "  Number of iterations: " << N_iter << "\n";
        std::cout<< "  Iteration matrix spectral radius: " << rho << "\n";
    }
}
