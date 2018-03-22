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
#include "Objects.h"


//==========================================================================
// Parser tools
//==========================================================================

// Parsing space(s) delimited strings into vector
template<class T>
std::vector<T> parse_vector(std::string const& pointLine)
{
    std::istringstream iss(pointLine);

    return std::vector<T>{ std::istream_iterator<T>(iss), 
                           std::istream_iterator<T>()
                         };
}

// Find object pointer by id number
template< typename T >
std::shared_ptr<T> find_by_id( const std::vector<std::shared_ptr<T>>& vec,
                               const int id )
{
    for ( auto& v : vec ){
	if ( v->id() == id ) { return v; }
    }
    return nullptr;
}


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

// The GLR algorithm to compute quadrature sets
void legendre_compute_glr ( int n, double x[], double w[] );



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
            std::vector<double> psi_b( N/2, 0.0 );
            for( int n = 0; n < N; n++ ){
                if( mu[n] - w[n]/2 <= bc_mu && bc_mu <= mu[n] + w[n]/2 ){
                    const double val = bc_mu * magnitude / mu[n] / w[n];
                    if( bc_side == "left" ){ psi_b[n-N/2] = val; }
                    if( bc_side == "right" ){ psi_b[n] = val; }
                    break;
                }
            }
            BC_set = std::make_shared<BCMonoDirectional>(psi_b);
        }
        if( bc_side == "left" )  { BC_left  = BC_set; }
        if( bc_side == "right" ) { BC_right = BC_set; }
    }


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
            S[j] = 0.5 * ( mesh[j]->SigmaS() * phi[j] + mesh[j]->Q() );
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
                psi[idx] = ( (mu[n] - 0.5 * mesh[j]->SigmaT() * mesh[j]->dz() ) 
                             * psi[idx] + S[j] * mesh[j]->dz() )
                           / ( mu[n] + 0.5 * mesh[j]->SigmaT() * mesh[j]->dz());
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
                psi[n] = ( ( -mu[n] - 0.5 * mesh[j]->SigmaT() * mesh[j]->dz() ) 
                           * psi[n] + S[j] * mesh[j]->dz() )
                         / ( -mu[n] + 0.5 * mesh[j]->SigmaT() * mesh[j]->dz() );
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

    } while ( error > ( 1.0 - rho.back() ) * epsilon );

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
    
    // Spectral radius estimates
    rho.erase(rho.begin());
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
