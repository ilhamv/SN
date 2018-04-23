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
                       const std::string accelerator_type,
                       const double beta )
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
    std::vector<double> phi_old(J);
    double              error;         // Maximum relative error
    double rho_num, rho_denom = 1.0;
    
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
                                                         BC_right, beta );
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
        BC_left->set_boundary( psi[0], 0.5*N, N );
        
        // Space sweep
        for( int j = 0; j < J; j++ ){

            // Set isotropic source
            S = 0.5 * ( mesh[j]->SigmaS() * phi[j] + mesh[j]->Q() );

            // Direction sweep
            for( int n = 0.5*N; n < N; n++ ){
                psi[j+1][n] = ( ( mu[n] - ( 1.0 - mesh[j]->alpha(n) ) * 0.5 
                                          * mesh[j]->tau() ) * psi[j][n]
                                + S * mesh[j]->dz() )
                              / ( mu[n] + ( 1.0 + mesh[j]->alpha(n) ) * 0.5 
                                          * mesh[j]->tau() );
            }
        }
        
        //======================================================================
        // Backward sweep
        //======================================================================

        // Set BC
        BC_right->set_boundary( psi[J], 0, 0.5*N );

        // Space sweep
        for( int j = J-1; j >= 0; j-- ){

            // Set isotropic source
            S = 0.5 * ( mesh[j]->SigmaS() * phi[j] + mesh[j]->Q() );

            // Direction sweep
            for( int n = 0; n < 0.5*N; n++ ){
                psi[j][n]   = ( ( -mu[n] - ( 1.0 + mesh[j]->alpha(n) ) * 0.5 
                                           * mesh[j]->tau() ) * psi[j+1][n]
                                + S * mesh[j]->dz() )
                              / ( -mu[n] + ( 1.0 - mesh[j]->alpha(n) ) * 0.5 
                                         * mesh[j]->tau() );
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
            double val = phi_old[j] / ( std::abs(phi[j]) + epsilon*epsilon );
// (commented out)( std::abs(phi[j]) + std::numeric_limits<double>::epsilon() );
            if( val > error ){ error = val; }
        }
        rho_num = norm_2(phi_old);
        rho.push_back( rho_num/rho_denom );
        rho_denom = rho_num;

    } while ( error > epsilon );
// (commented out)} while ( error > ( 1.0 - rho.back() ) * epsilon );
        
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

    // Space method, and set region's alpha
    std::cout<<"Space discretization: DD\n";

    // Source iteration tools
    const int N = mu.size();
    const int J = mesh.size();
    const int K = time.size() - 1;
    std::vector<std::vector<double>> psi;     // Cell-edge angular flux
    std::vector<std::vector<double>> psi_prv; // Previous
    double rho;                               // Spectral radius estimate
    double S;                                 // Isotropic source
    std::vector<double> phi_old(J);
    double error;                             // Maximum relative error
    double rho_num, rho_denom = 1.0;

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

    // Initialize phi
    phi.resize( K+1, std::vector<double>(J,0.0) );
    
    // Initial cell-edge angular flux
    psi = psi_initial;

    // Initial cell-average scalar flux
    for( int j = 0; j < J; j++ ){
        for( int n = 0; n < N; n++ ){
            phi[0][j] += ( psi[j][n] + psi[j+1][n] ) * w[n];
        }
        phi[0][j] *= 0.5;
    }


    //=========================================================================
    // March in time
    //=========================================================================
    
    for( int k = 1; k < K+1; k++ ){
        // Iteration counter
        int N_iter = 0;

        // Set first guess and previous time
        psi_prv = psi;
        phi[k] = phi[k-1];


        //======================================================================
        // Source Iteration - Implicit time
        //======================================================================
        
        do{
            //==================================================================
            // Forward sweep
            //==================================================================

            // Set BC
            BC_left->set_boundary( psi[0], 0.5*N, N );
            
            // Space sweep
            for( int j = 0; j < J; j++ ){
                
                // Isotropic source
                S = 0.5 * ( mesh[j]->SigmaS() * phi[k][j] + mesh[j]->Q() );
                
                // Direction sweep
                for( int n = 0.5*N; n < N; n++ ){
                    psi[j+1][n] = ( ( mu[n] - 0.5 * mesh[j]->tau() ) * psi[j][n]
                                    + ( S + aug * 0.5 
                                        * ( psi_prv[j][n] + psi_prv[j+1][n] ) )
                                    * mesh[j]->dz() )
                               / ( mu[n] + 0.5 * mesh[j]->tau() );
                }
            }
            
            //==================================================================
            // Backward sweep
            //==================================================================

            // Set BC
            BC_right->set_boundary( psi[J], 0, 0.5*N );

            // Space sweep
            for( int j = J-1; j >= 0; j-- ){
                
                // Isotropic source
                S = 0.5 * ( mesh[j]->SigmaS() * phi[k][j] + mesh[j]->Q() );
                
                for( int n = 0; n < 0.5*N; n++ ){
                    psi[j][n] = ( ( -mu[n] - 0.5 * mesh[j]->tau() )* psi[j+1][n]
                                  + ( S + aug * 0.5 
                                      * ( psi_prv[j][n] + psi_prv[j+1][n] ) )
                                  * mesh[j]->dz() )
                             / ( -mu[n] + 0.5 * mesh[j]->tau() );
                }
            }

            //==================================================================
            // Update flux
            //==================================================================
            
            // Store old phi
            phi_old = phi[k]; 

            // Update phi
            for( int j = 0; j < J; j++ ){
                phi[k][j] = 0.0;
                for( int n = 0; n < N; n++ ){
                    phi[k][j] += ( psi[j][n] + psi[j+1][n] ) * w[n];
                }
                phi[k][j] *= 0.5;
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
    
    // Space method, and set region's alpha
    std::cout<<"Space discretization: DD\n";

    // Source iteration tools
    const int N = mu.size();
    const int J = mesh.size();
    const int K = time.size() - 1;
    std::vector<std::vector<double>> psi;     // Cell-edge angular flux
    std::vector<std::vector<double>> psi_prv; // Previous
    double rho;                               // Spectral radius estimate
    double              psi_b;                // Boundary angular flux (t-avg)
    std::vector<double> phi_old(J);
    std::vector<double> Jc(J+1);              // Cell-edge current
    double              S,S0,S1;              // Isotropic sources
    double error;                             // Maximum relative error
    double rho_num, rho_denom = 1.0;

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
    std::vector<double> Ab(N), Bb(N), lambda(J);

    eta.resize( J, std::vector<double>(N,0.0) );
    xi.resize( J, std::vector<double>(N,0.0) );
    alpha.resize( J, std::vector<double>(N,0.0) );
    beta.resize( J, std::vector<double>(N,0.0) );
    A.resize( J, std::vector<double>(N,0.0) );
    B.resize( J, std::vector<double>(N,0.0) );
    C.resize( J, std::vector<double>(N,0.0) );

    for( int j = 0; j < J; j++ ){
        lambda[j] = mesh[j]->dz() * ( 0.5 * mesh[j]->SigmaA() + gamma );
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

    // Initialize phi
    phi.resize( K+1, std::vector<double>(J,0.0) );
    
    // Initial cell-edge angular flux
    psi = psi_initial;

    // Initial cell-average scalar flux
    for( int j = 0; j < J; j++ ){
        for( int n = 0; n < N; n++ ){
            phi[0][j] += ( psi[j][n] + psi[j+1][n] ) * w[n];
        }
        phi[0][j] *= 0.5;
    }
    
    // Initial cell-edge current
    for( int j = 0; j < J+1; j++ ){
        Jc[j] = 0.0;
        for( int n = 0; n < N; n++ ){
            Jc[j] += mu[n] * psi[j][n]* w[n];
        }
    }
    

    //=========================================================================
    // March in time
    //=========================================================================
    
    for( int k = 1; k < K+1; k++ ){
        // Iteration counter
        int N_iter = 0;

        // Set first guess and previous time
        psi_prv = psi;
        phi[k] = phi[k-1];

        //======================================================================
        // Source Iteration - Implicit time dependent
        //======================================================================
        
        do{

            //==================================================================
            // Begin forward sweep
            //==================================================================
            
            // Isotropic source
            S1 = mesh[0]->SigmaS() * phi[k][0] + mesh[0]->Q();
            S = 0.5 * mesh[0]->SigmaS() * ( Jc[1] - Jc[0] 
                                              + mesh[0]->tau() * phi[k][0] )
                + lambda[0] * S1;

            // Set BC
            BC_left->set_boundary( psi[0], 0.5*N, N );
            
            for( int n = 0.5*N; n < N; n++ ){
                
                // Set time-average BC
                if( BC_left->type() == "Reflective" ){
                    psi_b = 1.0 / ( 4.0 * mu[n] * gamma )
                        * (-( xi[0][n] * beta[0][n] + mesh[0]->dz() * gamma_sq )
                             * psi[1][N-1-n] 
                           -( xi[0][n] * alpha[0][n]+ mesh[0]->dz() * gamma_sq )
                             * psi[0][N-1-n] 
                           + S - mu[n] * S1 + mesh[0]->dz() * gamma_sq 
                           * ( psi_prv[0][N-1-n] + psi_prv[1][N-1-n] ) );
                } else{ psi_b = psi[0][n]; }

                psi[1][n] = ( S + mu[n] * S1  - Bb[n] * psi[0][n] 
                                + mesh[0]->dz() * gamma_sq 
                                  * ( psi_prv[0][n] + psi_prv[1][n] )
                                + 4.0 * mu[n] * gamma * psi_b ) / Ab[n];
            }

            //==================================================================
            // Continue forward sweep
            //==================================================================

            // Space sweep
            for( int j = 0; j < J-1; j++ ){

                // Isotropic source
                S0 = S1;
                S1 = mesh[j+1]->SigmaS() * phi[k][j+1] + mesh[j+1]->Q();
                S = 0.5* ( mesh[j+1]->SigmaS()
                           * (Jc[j+2] - Jc[j+1]+ mesh[j+1]->tau() * phi[k][j+1])
                           + mesh[j]->SigmaS()
                             * ( Jc[j+1] - Jc[j] + mesh[j]->tau() * phi[k][j] ))
                    + lambda[j+1] * S1 + lambda[j] * S0;

                // Direction sweep
                for( int n = 0.5*N; n < N; n++ ){
                    psi[j+2][n] = ( S + mu[n] * ( S1 - S0 )
                                    - B[j][n] * psi[j+1][n] 
                                    - C[j][n] * psi[j][n] 
                                    + gamma_sq * ( mesh[j+1]->dz()
                                       * ( psi_prv[j+1][n] + psi_prv[j+2][n] ) 
                                      + mesh[j]->dz()
                                       * ( psi_prv[j][n] + psi_prv[j+1][n] ) ) )
                                      / A[j][n];
                }
            }
            
            //==================================================================
            // Begin backward sweep
            //==================================================================
            
            // Isotropic source
            S0 = mesh[J-1]->SigmaS() * phi[k][J-1] + mesh[J-1]->Q();
            S = 0.5 * mesh[J-1]->SigmaS() * ( Jc[J] - Jc[J-1] 
                                              + mesh[J-1]->tau() * phi[k][J-1] )
                + lambda[J-1] * S0;

            // Set BC
            BC_right->set_boundary( psi[J], 0, 0.5*N );
            
            for( int n = 0; n < 0.5*N; n++ ){
                
                // Set time-average BC
                if( BC_right->type() == "Reflective" ){
                    psi_b = - 1.0 / ( 4.0 * mu[n] * gamma )
                        * ( -( eta[J-1][n] * beta[J-1][n] + mesh[J-1]->dz() 
                                                            * gamma_sq ) 
                             * psi[J][N-1-n]
                            -( eta[J-1][n] * alpha[J-1][n] + mesh[J-1]->dz() 
                                                             * gamma_sq ) 
                             * psi[J-1][N-1-n]
                           + S + mu[n] * S0
                           + mesh[J-1]->dz() * gamma_sq 
                             * ( psi_prv[J-1][N-1-n] + psi_prv[J][N-1-n] ) );
                } else{ psi_b = psi[J][n]; }

                psi[J-1][n] = ( S - mu[n] * S0 - Ab[n] * psi[J][n]
                                    + mesh[J-1]->dz() * gamma_sq 
                                      * ( psi_prv[J-1][n] + psi_prv[J][n] )
                                    - 4.0 * mu[n] * gamma * psi_b ) / Bb[n];
            }

            //==================================================================
            // Continue backward sweep
            //==================================================================

            for( int j = J-2; j >= 0; j-- ){

                // Isotropic source
                S1 = S0;
                S0 = mesh[j]->SigmaS() * phi[k][j] + mesh[j]->Q();
                S = 0.5* ( mesh[j+1]->SigmaS()
                           * (Jc[j+2] - Jc[j+1]+ mesh[j+1]->tau() * phi[k][j+1])
                           + mesh[j]->SigmaS()
                             * ( Jc[j+1] - Jc[j] + mesh[j]->tau() * phi[k][j] ))
                    + lambda[j+1] * S1 + lambda[j] * S0;

                for( int n = 0; n < 0.5*N; n++ ){
                    psi[j][n] = ( S + mu[n] * ( S1 - S0 )
                                      - A[j][n] * psi[j+2][n]
                                      - B[j][n] * psi[j+1][n] 
                                      + gamma_sq * ( mesh[j+1]->dz() 
                                         * ( psi_prv[j+1][n] + psi_prv[j+2][n] )
                                        + mesh[j]->dz() 
                                       * ( psi_prv[j][n] + psi_prv[j+1][n] ) ) )
                                      / C[j][n];
                }
            }

            //==================================================================
            // Update cell-average flux and cell-edge current
            //==================================================================
            
            // Store old phi
            phi_old = phi[k]; 

            // Update phi
            for( int j = 0; j < J; j++ ){
                phi[k][j] = 0.0;
                for( int n = 0; n < N; n++ ){
                    phi[k][j] += ( psi[j][n] + psi[j+1][n] ) * w[n];
                }
                phi[k][j] *= 0.5;
            }

            for( int j = 0; j < J+1; j++ ){
                Jc[j] = 0;
                for( int n = 0; n < N; n++ ){
                    Jc[j] += mu[n] * psi[j][n]* w[n];
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

        } while ( error > ( 1.0 - rho) * epsilon );
        
        // Some outputs
        std::cout<< "Report for k = " << k << " ("<< time[k] << " s)\n";
        std::cout<< "  Number of iterations: " << N_iter << "\n";
        std::cout<< "  Iteration matrix spectral radius: " << rho << "\n";
    }
}
