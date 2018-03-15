#include <iostream>
#include <cstring> 
#include <iterator>
#include <vector>
#include <sstream>

#include "pugixml.hpp"

template<class T>
std::vector<T> parse_vector(std::string const& pointLine)
{
    std::istringstream iss(pointLine);

    return std::vector<T>{ std::istream_iterator<T>(iss), 
                           std::istream_iterator<T>()
                         };
}

int main( int argc, char* argv[] )
{
    //==========================================================================
    // Parse input
    //==========================================================================

    // I/O directory
    if ( argc == 1 ){ 
        std::cout<< "[INPUT ERROR] Please provide I/O directory...\n";
        std::exit(EXIT_FAILURE);
    }
    const std::string io_dir = std::string(argv[1]) + "/";

    // XML input file
    std::string input_name = io_dir + "input.xml";
    pugi::xml_document input_file;
    input_file.load_file(input_name.c_str());

    // Quadrature and convergence criterion
    const int          N = std::stoi( input_file.child_value("quadrature") );
    const double epsilon = std::stod( input_file.child_value("convergence") );

    // Region edges and number of meshes
    const std::vector<double> R_edge  = parse_vector<double>
                                        ( input_file.child_value("edge") );
    const std::vector<int>    R_mesh = parse_vector<int>
                                       ( input_file.child_value("mesh") );
    const int N_region = R_mesh.size();

    // Region properties
    const std::vector<double> R_SigmaT = parse_vector<double>
                                         ( input_file.child_value("total") );
    const std::vector<double> R_SigmaS = parse_vector<double>
                                         ( input_file.child_value("scatter") );
    const std::vector<double> R_Q      = parse_vector<double>
                                         ( input_file.child_value("source") );

    // Boundary conditions
    const std::vector<int> bc = parse_vector<int>
                                ( input_file.child_value("bc") );


    //==========================================================================
    // Construct the meshes
    //==========================================================================

    int idx; // Index helper

    // # of meshes
    int N_mesh = 0;
    for( int i = 0; i < N_region; i++ ){ N_mesh += R_mesh[i]; }

    // Mesh size
    idx = 0;
    std::vector<double> dz(N_mesh);
    for( int i = 0; i < N_region; i++ ){
        for( int j = 0; j < R_mesh[i]; j++ ){
            dz[idx] = ( R_edge[i+1] - R_edge[i] ) / R_mesh[i];
            idx++;
        }
    }

    // Mesh properties
    std::vector<double> SigmaT(N_mesh);
    std::vector<double> SigmaS(N_mesh);
    std::vector<double> Q(N_mesh);
    idx = 0;
    for( int i = 0; i < N_region; i++ ){
        for( int j = 0; j < R_mesh[i]; j++ ){
            SigmaT[idx] = R_SigmaT[i];
            SigmaS[idx] = R_SigmaS[i];
            Q[idx]      = R_Q[i];
            idx++;
        }
    }


    return 0;
}
