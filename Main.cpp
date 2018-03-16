#include <iostream>
#include <cstring> 
#include <iterator>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>

#include "pugixml.hpp"
#include "H5Cpp.h"


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
        val = std::abs( phi[j] - phi_old[j] ) / phi[j];
        if( val > emax ){ emax = val; }
    }   
    return emax;
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
    const std::vector<std::string> bc = parse_vector<std::string>
                                        ( input_file.child_value("bc") );
    std::vector<bool> reflect(2,false);
    if( bc[0] == "reflective" ) { reflect[0] = true; }
    if( bc[1] == "reflective" ) { reflect[1] = true; }


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

    // Angular flux at boundary (default: vacuum)
    std::vector<double> psi_left( N/2, 0.0 );
    std::vector<double> psi_right( N/2, 0.0 );

    if( bc[0] == "vacuum" ) { ; } // Skip
    if( bc[1] == "vacuum" ) { ; } // Skip


    //==========================================================================
    // The Source Iteration
    //==========================================================================

    std::vector<double> psi( N/2 );    // In/out angular flux
    double psi_avg;                    // Average angular flux
    std::vector<double> phi( J, 0.0 ); // Scalar flux, zero first guess
    std::vector<double> phi_old(J);
    std::vector<double> S(J);          // RHS source

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

        // Forward sweep
        psi = psi_left;
        for( int j = 0; j < J; j++ ){
            idx = 0; // Index for psi, as its size is N/2
            for( int n = 0.5*N; n < N; n++ ){
                psi_avg = psi[idx];
                psi[idx] = ( ( std::abs(mu[n]) - 0.5 * SigmaT[j] * dz[j] ) * psi[idx]
                             + S[j] * dz[j] )
                           / ( std::abs(mu[n]) + 0.5 * SigmaT[j] * dz[j] );
                psi_avg = (psi_avg + psi[idx]) * w[n];
                phi[j] += psi_avg;
                idx++;
            }
        }
        
        // Backward sweep
        //psi = psi_right;
        std::reverse(psi.begin(), psi.end());
        for( int j = J-1; j >= 0; j-- ){
            for( int n = 0; n < 0.5*N; n++ ){
                psi_avg = psi[n];
                psi[n] = ( ( std::abs(mu[n]) - 0.5 * SigmaT[j] * dz[j] ) * psi[n]
                           + S[j] * dz[j] )
                         / ( std::abs(mu[n]) + 0.5 * SigmaT[j] * dz[j] );
                psi_avg = (psi_avg + psi[n]) * w[n];
                phi[j] += psi_avg;
            }
            phi[j] *= 0.5;
        }
    } while( maximum_relative_error( phi, phi_old ) > epsilon );


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

    dataset = output.createDataSet( "N", type_int, space_scalar );
    dataset.write(&N, type_int);
    dataset = output.createDataSet( "epsilon", type_double, space_scalar );
    dataset.write(&epsilon, type_double);
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
    dataset.write(bc[0], type_string);
    dataset = output.createDataSet( "bc_right", type_string, space_scalar);
    dataset.write(bc[1], type_string);

    dims[0] = N;
    space_vector = H5::DataSpace(1,dims);
    dataset = output.createDataSet( "mu_n", type_double, space_vector);
    dataset.write(mu.data(), type_double);
    dataset = output.createDataSet( "w_n", type_double, space_vector);
    dataset.write(w.data(), type_double);

    dims[0] = J;
    space_vector = H5::DataSpace(1,dims);
    dataset = output.createDataSet( "dz", type_double, space_vector);
    dataset.write(dz.data(), type_double);
    dataset = output.createDataSet( "tot", type_double, space_vector);
    dataset.write(SigmaT.data(), type_double);
    dataset = output.createDataSet( "scat", type_double, space_vector);
    dataset.write(SigmaS.data(), type_double);
    dataset = output.createDataSet( "source", type_double, space_vector);
    dataset.write(Q.data(), type_double);

    // Scalar flux
    dims[0] = J;
    space_vector = H5::DataSpace(1,dims);
    dataset = output.createDataSet( "scalar_flux", type_double, space_vector);
    dataset.write(phi.data(), type_double);
    dataset = output.createDataSet( "z", type_double, space_vector);
    dataset.write(z.data(), type_double);

    return 0;
}
