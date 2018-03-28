#include <iostream>
#include <cstring> 
#include <vector>
#include <limits> 
#include <memory>

#include "pugixml.hpp"
#include "H5Cpp.h"

#include "Objects.h"
#include "Algorithm.h"
#include "Solver.h"


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
    // Method descriptions
    //==========================================================================

    // Quadrature
    const int N = std::stoi( input_file.child_value("N") );
    
    // Convergence criterion
    const double epsilon = std::stod( input_file.child_value("epsilon") );

    // Spatial discretization method
    const std::string method = input_file.child_value("method");


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

    // Set up region properties
    for( int i = 0; i < N_region; i++ ){
        const double r_dz = r_space[i] / r_mesh[i];
        const std::shared_ptr<Material> r_M = 
                                          find_by_id( material, r_material[i] );
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
    // Quadrature sets
    //==========================================================================

    std::vector<double> mu(N);
    std::vector<double> w(N);
    
    // Calling the GLR algorithm
    legendre_compute_glr( N, &mu[0], &w[0]);
   

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
                    if( bc_side == "left" ){ b_idx = n-N/2; }
                    if( bc_side == "right" ){ b_idx = n; }
                    break;
                }
            }
            BC_set = std::make_shared<BCMonoDirectional>(b_val, b_idx);
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
    std::vector<double> time = {0.0};
    double dt;
    double speed;
    int K;

    if( input_file.child("TD") ){
        TD = true;
        pugi::xml_node input_TD = input_file.child("TD");

        // Initial condition
        std::string ic_type = input_TD.child("IC").attribute("type").value();
        if( ic_type == "zero" ){
            psi_initial.resize(J+1, std::vector<double>(N,0.0));
        }
        
        // Time step
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
    // Steady State Solver
    //==========================================================================

    std::vector<double> phi; // Cell-average scalar flux, with zero first guess
    std::vector<double> rho; // Spectral radius
    int N_iter;

    if( !TD ){
        // Initialize phi
        phi.resize( J, 0.0 );

        // Solve!
        source_iteration( epsilon, mesh, mu, w, BC_left, BC_right, phi, rho );

        // Some outputs
        N_iter = rho.size();
        std::cout<< "Done!\n";
        std::cout<< "Number of iterations: " << N_iter << "\n";
        std::cout<< "Iteration matrix spectral radius: " << rho.back() << "\n";
    }
    

    //==========================================================================
    // Time Dependent Solver
    //==========================================================================

    std::vector<std::vector<double>> phi_t; // Cell-average scalar flux,
                                            //   at each time step

    if( TD ){
        // Time augment
        const double aug = 1.0 / speed / dt;

        // Initialize cell-average scalar flux at each time k
        phi_t.resize( K+1, std::vector<double>(J,0.0) );

        // Solve!
        source_iteration_TD( epsilon, mesh, region, mu, w, BC_left, BC_right,
                             speed, dt, K, psi_initial, phi_t );    
    }
    

    //==========================================================================
    // HDF5 output
    //==========================================================================
   
    H5std_string FILE_NAME(io_dir + "output.h5");
    H5::H5File output(FILE_NAME, H5F_ACC_TRUNC);
    H5::DataSet dataset;
    H5::Group group;
    H5::DataSpace space_scalar(H5S_SCALAR);
    H5::DataSpace space_vector;;
    H5::DataType type_int    = H5::PredType::NATIVE_INT;
    H5::DataType type_double = H5::PredType::NATIVE_DOUBLE;
    H5::StrType  type_string(0, H5T_VARIABLE);
    hsize_t dims[1]; 

    // Problem summary
    dataset = output.createDataSet( "N", type_int, space_scalar );
    dataset.write(&N, type_int);
    dataset = output.createDataSet( "epsilon", type_double, space_scalar );
    dataset.write(&epsilon, type_double);
    if(TD){
        dataset = output.createDataSet( "N_iter", type_double, space_scalar );
        dataset.write(&N_iter, type_double);
    }
    // Material, Region, Mesh
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
    }
    
    //==========================================================================
    // Time dependent outputs
    //==========================================================================

    if( TD ){
        phi.resize((K+1)*J);
        idx = 0;
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
        dataset = output.createDataSet( "scalar_flux", type_double, 
                                        space_vector);
        dataset.write(phi.data(), type_double);
        dataset = output.createDataSet( "z", type_double, space_vector);
        dataset.write(z.data(), type_double);
    }
    
    return 0;
}
