#include <iostream>
#include <cstring> 
#include <iterator>
#include <vector>
#include <sstream>
#include <cmath>
#include <limits> 
#include <memory>

#include "pugixml.hpp"
#include "H5Cpp.h"
#include "BC.h"


// The GLR algorithm to compute quadrature sets
void legendre_compute_glr ( int n, double x[], double w[] );

// Parsing space(s) delimited strings into vector
template<class T>
std::vector<T> parse_vector(std::string const& pointLine)
{
    std::istringstream iss(pointLine);

    return std::vector<T>{ std::istream_iterator<T>(iss), 
                           std::istream_iterator<T>()
                         };
}

// Maximum relative error of two vectors
double maximum_relative_error( const std::vector<double>& phi,
                               const std::vector<double>& phi_old )
{
    double emax = 0; double val;
    for( int j = 0; j < phi.size(); j++ ){
        val = std::abs( phi[j] - phi_old[j] ) 
              / ( std::abs(phi[j]) + std::numeric_limits<double>::epsilon() );
        if( val > emax ){ emax = val; }
    }   
    return emax;
}

double norm_2( const std::vector<double> v )
{
    double norm = 0.0;
    for( int i = 0; i < v.size(); i++ ){
        norm += v[i]*v[i];
    }
    norm = std::sqrt(norm);
    return norm;
}

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

    // Quadrature and convergence criterion
    const int          N = std::stoi( input_file.child_value("N") );
    const double epsilon = std::stod( input_file.child_value("epsilon") );

    // Region edges and number of meshes
    const std::vector<double> R_space = parse_vector<double>
                                        ( input_file.child_value("space") );
    const std::vector<int>    R_mesh  = parse_vector<int>
                                        ( input_file.child_value("mesh") );
    const int N_region = R_mesh.size();

    // Region properties
    const std::vector<double> R_SigmaT = parse_vector<double>
                                         ( input_file.child_value("SigmaT") );
    const std::vector<double> R_SigmaS = parse_vector<double>
                                         ( input_file.child_value("SigmaS") );
    const std::vector<double> R_Q      = parse_vector<double>
                                         ( input_file.child_value("Q") );

    // Boundary conditions
    std::shared_ptr<BC> BC_left, BC_right;
    
    for( auto bc : input_file.child("bc").children() ){
        std::shared_ptr<BC> BC_set;
        const std::string bc_type = bc.attribute("type").value();
        if( bc_type == "vacuum" ){
            BC_set = std::make_shared<BCVacuum>();
        } else if( bc_type == "reflective" ){
            BC_set = std::make_shared<BCReflective>();
        }
        if( std::string(bc.name()) == "left" )  { BC_left  = BC_set; }
        if( std::string(bc.name()) == "right" ) { BC_right = BC_set; }
    }


    //==========================================================================
    // The quadrature sets
    //==========================================================================

    std::vector<double> mu(N);
    std::vector<double> w(N);
    
    // Calling the GLR algorithm
    legendre_compute_glr( N, &mu[0], &w[0]);


    //==========================================================================
    // The meshes
    //==========================================================================

    int idx; // Index helper

    // # of meshes
    int J = 0;
    for( int i = 0; i < N_region; i++ ){ J += R_mesh[i]; }

    // Mesh size
    idx = 0;
    std::vector<double> dz(J);
    for( int i = 0; i < N_region; i++ ){
        for( int j = 0; j < R_mesh[i]; j++ ){
            dz[idx] = R_space[i] / R_mesh[i];
            idx++;
        }
    }

    // Center points
    std::vector<double> z(J);
    z[0] = 0.5 * dz[0];
    for( int j = 1; j < J; j++ ){
        z[j] = z[j-1] + 0.5 * ( dz[j-1] + dz[j] );
    }

    // Mesh properties
    std::vector<double> SigmaT(J);
    std::vector<double> SigmaS(J);
    std::vector<double> Q(J);
    idx = 0;
    for( int i = 0; i < N_region; i++ ){
        for( int j = 0; j < R_mesh[i]; j++ ){
            SigmaT[idx] = R_SigmaT[i];
            SigmaS[idx] = R_SigmaS[i];
            Q[idx]      = R_Q[i];
            idx++;
        }
    }

    
    //==========================================================================
    // Set angular flux B.C.
    //==========================================================================
/*
    // Angular flux at boundary (default: vacuum)
    std::vector<double> psi_left( N/2, 0.0 );
    std::vector<double> psi_right( N/2, 0.0 );

    // Reflective boundary is treated with toggle "reflect"

    // Isotropic
    for( int i = 0; i < 2; i++ ){ 
        if( bc[i] == "isotropic" ){
            for( int n = 0; n < N/2; n++ ){
                if( i == 0 ) { psi_left[n] = 1.0; }
                if( i == 1 ) { psi_right[n] = 1.0; }
            }
        }
    }

    // Monodirectional
    for( int i = 0; i < 2; i++ ){ 
        if( bc[i] != "vacuum" && bc[i] != "reflective" && bc[i] != "isotropic" )
        {
            const double psi_b = std::stod( bc[i] );
            for( int n = 0; n < N; n++ ){
                if( mu[n] - w[n]/2 <= psi_b && psi_b <= mu[n] + w[n]/2 ){
                    if( i == 0 ){ psi_left[n-N/2] = 1.0/w[n]; break;}
                    if( i == 1 ){ psi_right[n] = 1.0/w[n]; break;}
                }
            }
        }
    }
*/

    //==========================================================================
    // The Source Iteration
    //==========================================================================

    std::vector<double> psi( N/2 );    // In/out angular flux
    double psi_avg;                    // Average angular flux
    std::vector<double> phi( J, 0.0 ); // Scalar flux, zero first guess
    std::vector<double> phi_old(J);
    std::vector<double> S(J);          // RHS source
    double error;                      // Maximum relative error
    std::vector<double> rho;           // Spectral radius
    double rho_num, rho_denom = 1.0;

    // Start iteration
    do{

        // RHS source (note: isotropic)
        for( int j = 0; j < J; j++ ){
            S[j] = 0.5 * ( SigmaS[j] * phi[j] + Q[j] );
        }

        // Old phi
        phi_old = phi; 
        // Reset phi
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
                psi[idx] = ( (mu[n] - 0.5 * SigmaT[j] * dz[j] ) * psi[idx]
                             + S[j] * dz[j] )
                           / ( mu[n] + 0.5 * SigmaT[j] * dz[j] );
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
                psi[n] = ( ( -mu[n] - 0.5 * SigmaT[j] * dz[j] ) * psi[n]
                           + S[j] * dz[j] )
                         / ( -mu[n] + 0.5 * SigmaT[j] * dz[j] );
                psi_avg = (psi_avg + psi[n]) * w[n];
                phi[j] += psi_avg;
            }
            phi[j] *= 0.5;
        }
        
        //======================================================================
        // Relative error and spectral radius estimate
        //======================================================================
       
        error = 0.0;
        for( int j = 0; j < phi.size(); j++ ){
            // Now phi_old holds the absolute difference between iterates
            phi_old[j] = std::abs( phi[j] - phi_old[j] );
            double val = phi_old[j] / 
                  ( std::abs(phi[j]) + std::numeric_limits<double>::epsilon() );
            if( val > error ){ error = val; }
        }

        rho_num = norm_2(phi_old);
        rho.push_back( rho_num/rho_denom );
        rho_denom = rho_num;

    } while ( error > epsilon );

    // Some outputs
    const int N_iter = rho.size();
    std::cout<< "Done!\n";
    std::cout<< "Number of iterations: " << N_iter << "\n";
    std::cout<< "Iteration matrix spectral radius: " << rho.back() << "\n";
    

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
    dataset = output.createDataSet( "N_iter", type_double, space_scalar );
    dataset.write(&N_iter, type_double);
    dims[0] = N_region;
    space_vector = H5::DataSpace(1,dims);
    dataset = output.createDataSet( "space", type_double, space_vector);
    dataset.write(R_space.data(), type_double);
    dataset = output.createDataSet( "mesh", type_int, space_vector);
    dataset.write(R_mesh.data(), type_int);
    dataset = output.createDataSet( "SigmaT", type_double, space_vector);
    dataset.write(R_SigmaT.data(), type_double);
    dataset = output.createDataSet( "SigmaS", type_double, space_vector);
    dataset.write(R_SigmaS.data(), type_double);
    dataset = output.createDataSet( "Q", type_double, space_vector);
    dataset.write(R_Q.data(), type_double);
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
    
    // Spectral radius estimates
    rho[0] = 1.0;
    dims[0] = rho.size();
    space_vector = H5::DataSpace(1,dims);
    dataset = output.createDataSet( "spectral_radius",type_double,space_vector);
    dataset.write(rho.data(), type_double);

    // Scalar flux solution
    dims[0] = J;
    space_vector = H5::DataSpace(1,dims);
    dataset = output.createDataSet( "scalar_flux", type_double, space_vector);
    dataset.write(phi.data(), type_double);
    dataset = output.createDataSet( "z", type_double, space_vector);
    dataset.write(z.data(), type_double);

    return 0;
}
