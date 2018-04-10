#include <iostream>
#include <cstring> 
#include <vector>
#include <limits> 
#include <memory>
#include <cmath>

#include "pugixml.hpp"
#include "H5Cpp.h"

#include "Objects.h"
#include "Algorithm.h"
#include "Solver.h"
#include "Accelerator.h"


int main( int argc, char* argv[] )
{
    //==========================================================================
    // Parse XML input
    //==========================================================================

    // I/O directory
    const std::string io_dir = std::string(argv[1]) + "/";

    // XML input file
    std::string input_name = io_dir + "input.xml";
    pugi::xml_document input_file;
    input_file.load_file(input_name.c_str());


    //==========================================================================
    // Method parameters
    //==========================================================================

    // Quadrature
    const int N = std::stoi( input_file.child_value("N") );
    
    // Convergence criterion
    const double epsilon = std::stod( input_file.child_value("epsilon") );
    
    // Accelerator
    const std::string accelerator_type = input_file.child("Accelerator").
                                                    attribute("type").value();


    //==========================================================================
    // Quadrature sets
    //==========================================================================

    std::vector<double> mu(N);
    std::vector<double> w(N);
    
    // Calling the GLR algorithm
    legendre_compute_glr( N, &mu[0], &w[0]);
   

    //==========================================================================
    // Materials
    //==========================================================================
    
    std::vector<std::shared_ptr<Material>> material;

    for( auto m : input_file.child("materials").children("material") ){
        const int         m_id      = m.attribute("id").as_int();
        const std::string m_name    = m.attribute("name").value();
        const double      m_total   = m.child("total").
                                        attribute("xs").as_double();
        const double      m_scatter = m.child("scatter").
                                        attribute("xs").as_double();
        material.push_back( std::make_shared<Material>( m_id, m_name, 
                                                        m_total, m_scatter ) );
    }


    //==========================================================================
    // Regions
    //==========================================================================
    
    std::vector<std::shared_ptr<Region>> region;
    
    pugi::xml_node input_region = input_file.child("region");

    // Region space, number of meshes, material number, and source strength
    const std::vector<double> r_space = parse_vector<double>
                                        ( input_region.child_value("space") );
    const std::vector<int> r_mesh = parse_vector<int>
                                    ( input_region.child_value("mesh") );
    const std::vector<int> r_material = parse_vector<int>
                                        ( input_region.child_value("material"));
    const std::vector<double> r_Q = parse_vector<double>
                                    ( input_region.child_value("source") );
    const int N_region = r_space.size();

    // Space discretization method (Steady state only)
    std::string space_method = input_region.attribute("method").value();

    // Set up region properties
    for( int i = 0; i < N_region; i++ ){
        // Mesh size
        const double r_dz = r_space[i] / r_mesh[i];
        // Material
        const std::shared_ptr<Material> r_M = 
                                          find_by_id( material, r_material[i] );
        // Create region
        region.push_back( std::make_shared<Region>( r_M, r_dz, r_Q[i] ));
    }


    //==========================================================================
    // Meshes
    //==========================================================================

    std::vector<std::shared_ptr<Region>> mesh;

    int idx; // Index helper

    // # of meshes
    int J = 0;
    for( int i = 0; i < N_region; i++ ){ J += r_mesh[i]; }
    mesh.resize(J);

    // Point mesh to the corresponding region
    idx = 0;
    for( int i = 0; i < N_region; i++ ){
        for( int j = 0; j < r_mesh[i]; j++ ){
            mesh[idx] = region[i];
            idx++;
        }
    }

    // Center points
    std::vector<double> z(J);
    z[0] = 0.5 * mesh[0]->dz();
    for( int j = 1; j < J; j++ ){
        z[j] = z[j-1] + 0.5 * ( mesh[j-1]->dz() + mesh[j]->dz() );
    }


    //==========================================================================
    // Boundary conditions
    //==========================================================================

    std::shared_ptr<BC> BC_left, BC_right;
    
    for( auto bc : input_file.child("bc").children() ){
        std::shared_ptr<BC> BC_set;
        const std::string bc_type = bc.attribute("type").value();
        const std::string bc_side = bc.name();

        if( bc_type == "vacuum" ){
            BC_set = std::make_shared<BCVacuum>();
        } else if( bc_type == "reflective" ){
            BC_set = std::make_shared<BCReflective>();
        } else if( bc_type == "isotropic" ){
            const double magnitude = bc.attribute("magnitude").as_double();
            BC_set = std::make_shared<BCIsotropic>(magnitude);
        } else if( bc_type == "mono" ){
            const double magnitude = bc.attribute("magnitude").as_double();
            const double bc_mu     = bc.attribute("mu").as_double();
            double b_val;
            unsigned long long b_idx;
            for( int n = 0; n < N; n++ ){
                if( mu[n] - w[n]/2 <= bc_mu && bc_mu <= mu[n] + w[n]/2 ){
                    b_val = bc_mu * magnitude / mu[n] / w[n];
                    b_idx = n;
                    break;
                }
            }
            BC_set = std::make_shared<BCMonoDirectional>(b_val, b_idx);
        } else if( bc_type == "linear" ){
            std::vector<double> bc_param = parse_vector<double>( bc.
                                                attribute("param").value() );
            const double a = bc_param[0];
            const double b = bc_param[1];
            std::vector<double> psi_b;
            int idx; double val, mu1, mu2;
            if( bc_side == "left" )  { idx = N/2; }
            if( bc_side == "right" ) { idx = 0; }
            for( int n = idx; n<idx+N/2; n++ ){
                mu1 = -1.0;
                for( int m = 0; m < n; m++ ){ mu1 += w[m]; }
                mu2 = mu1 + w[n];
                std::cout<<n<<"  "<<mu1<<"  "<<mu2<<"\n";
                val = 1.0 / ( mu[n]*w[n] ) 
                      * ( 0.5 * a * mu2*mu2 + 1.0/3.0 * b * mu2*mu2*mu2
                          - ( 0.5 * a * mu1*mu1 + 1.0/3.0 * b * mu1*mu1*mu1 ) );
                psi_b.push_back(std::abs(val));
            }
            BC_set = std::make_shared<BCLinear>(psi_b);
        }
        if( bc_side == "left" )  { BC_left  = BC_set; }
        if( bc_side == "right" ) { BC_right = BC_set; }
    }


    //==========================================================================
    // Time dependent input
    //==========================================================================

    bool TD = false;
    std::vector<std::vector<double>> psi_initial; // Initial cell-edge angular 
                                                  //   flux [J][N]
    const std::string TD_method = input_file.child("TD").attribute("method")
                                                        .value();
    std::vector<double> time = {0.0};
    double dt;
    double speed;
    int K;

    if( input_file.child("TD") ){
        TD = true;
        pugi::xml_node input_TD = input_file.child("TD");

        // Initial condition
        std::string ic_type = input_TD.child("IC").attribute("type").value();
        std::vector<double> ic_param = parse_vector<double>( input_TD.
                                       child("IC").attribute("param").value() );
        const std::vector<double> r_space = parse_vector<double>
                                        ( input_region.child_value("space") );
        if( ic_type == "zero" ){
            psi_initial.resize(J+1, std::vector<double>(N,0.0));
        } else if( ic_type == "one" ){
            // alpha + beta mu^2
            double a,b,alpha,beta;
            std::vector<double> psi_ic(N);
            alpha = ic_param[0]; beta = ic_param[1];

            for( int n = 0; n < N; n++ ){
                a = -1.0;
                for( int m = 0; m < n; m++ ){ a += w[m]; }
                b = a + w[n];
                psi_ic[n] = 1.0 / w[n] * ( alpha * w[n] 
                                           + beta / 3.0 * ( b*b*b - a*a*a ) );
            }
            psi_initial.resize(J+1, psi_ic);
        }
        
        // Time step and time
        pugi::xml_node input_time = input_file.child("TD").child("time");
        K = input_time.attribute("step").as_int();
        dt = input_time.attribute("final").as_double() / K;
        for( int k = 0; k < K; k++ ){
            time.push_back(time.back() + dt);
        }

        // Speed
        speed = std::stod( input_TD.child_value("speed") );
    }


    //==========================================================================
    // Solve: Steady state
    //==========================================================================

    // Results
    std::vector<double> phi;              // Cell-average scalar flux
    std::vector<std::vector<double>> psi; // Cell-edge angular flux
    std::vector<double> rho;              // Spectral radius
    int                 N_iter;           // # of iterations

    if( !TD ){
        source_iteration( N_iter, epsilon, mesh, region, mu, w, BC_left, 
                          BC_right, phi, psi, rho, space_method, 
                          accelerator_type);
    }
    

    //==========================================================================
    // Solve: Time dependent
    //==========================================================================

    // Results
    std::vector<std::vector<double>> phi_t; // Cell-average scalar flux,
                                            //   at each time step

    if( TD ){
        if( TD_method == "implicit" ){
            source_iteration_TD_implicit( epsilon, mesh, material, region, mu, 
                                          w, BC_left, BC_right, speed, dt, 
                                          psi_initial, phi_t, accelerator_type, 
                                          time );        
        }
        if( TD_method == "MB" ){
            source_iteration_TD_MB( epsilon, mesh, material, region, mu, w, 
                                    BC_left, BC_right, speed, dt, psi_initial, 
                                    phi_t, accelerator_type, time );        
        }
    }
    

    //==========================================================================
    // HDF5 output
    //==========================================================================
   
    H5std_string FILE_NAME(io_dir + "output.h5");
    H5::H5File output(FILE_NAME, H5F_ACC_TRUNC);
    H5::DataSet dataset;
    H5::Group group;
    H5::DataSpace space_scalar(H5S_SCALAR);
    H5::DataSpace space_vector;
    H5::DataSpace space_2D;
    H5::DataType type_int    = H5::PredType::NATIVE_INT;
    H5::DataType type_double = H5::PredType::NATIVE_DOUBLE;
    H5::StrType  type_string(0, H5T_VARIABLE);
    hsize_t dims[1]; 
    hsize_t dims_2D[2]; 

    // Problem summary
    dataset = output.createDataSet( "N", type_int, space_scalar );
    dataset.write(&N, type_int);
    dataset = output.createDataSet( "epsilon", type_double, space_scalar );
    dataset.write(&epsilon, type_double);
    
    // BC
    dataset = output.createDataSet( "bc_left", type_string, space_scalar);
    dataset.write(BC_left->type(), type_string);
    dataset = output.createDataSet( "bc_right", type_string, space_scalar);
    dataset.write(BC_right->type(), type_string);
    // Quadrature sets
    dims[0] = N;
    space_vector = H5::DataSpace(1,dims);
    dataset = output.createDataSet( "mu_n", type_double, space_vector);
    dataset.write(mu.data(), type_double);
    dataset = output.createDataSet( "w_n", type_double, space_vector);
    dataset.write(w.data(), type_double);
    
    //==========================================================================
    // Steady state outputs
    //==========================================================================

    if( !TD ){
        // # of iterations
        dataset = output.createDataSet( "N_iter", type_int, space_scalar );
        dataset.write(&N_iter, type_int);

        // Spectral radius estimates
        rho.erase(rho.begin());
        dims[0] = rho.size();
        space_vector = H5::DataSpace(1,dims);
        dataset = output.createDataSet( "spectral_radius",type_double,
                                        space_vector);
        dataset.write(rho.data(), type_double);

        // Scalar flux solution
        dims[0] = J;
        space_vector = H5::DataSpace(1,dims);
        dataset = output.createDataSet( "scalar_flux", type_double, 
                                        space_vector);
        dataset.write(phi.data(), type_double);
        dataset = output.createDataSet( "z", type_double, space_vector);
        dataset.write(z.data(), type_double);
        
        // Angular flux solution
        dims_2D[0] = J+1;
        dims_2D[1] = N;
        space_2D = H5::DataSpace(2,dims_2D);
        dataset = output.createDataSet( "angular_flux", type_double, 
                                        space_2D);
        phi.resize(N*(J+1));
        for( int j = 0; j < J+1; j++ ){
            for( int n = 0; n < N; n++ ){
                phi[N*j+n] = psi[j][n];
            }
        }
        dataset.write(phi.data(), type_double);
    }
    
    //==========================================================================
    // Time dependent outputs
    //==========================================================================

    if( TD ){
        phi.resize((K+1)*J);
        for( int k = 0; k < K+1; k++ ){
            for( int j = 0; j < J; j++ ){
                phi[J*k+j] = phi_t[k][j];
            }
        }
        // Scalar flux
        hsize_t dimsM[2]; dimsM[0] = K+1; dimsM[1] = J;
        H5::DataSpace data_spaceM(2,dimsM);
        dataset = output.createDataSet( "scalar_flux_time", type_double, 
                                        data_spaceM);
        dataset.write(phi.data(), type_double);
        
        // z
        dims[0] = J;
        space_vector = H5::DataSpace(1,dims);
        dataset = output.createDataSet( "z", type_double, space_vector);
        dataset.write(z.data(), type_double);
        
        // time
        dims[0] = K+1;
        space_vector = H5::DataSpace(1,dims);
        dataset = output.createDataSet( "time", type_double, 
                                        space_vector);
        dataset.write(time.data(), type_double);
    }
    
    return 0;
}
