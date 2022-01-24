#include <iostream>
#include "cblas.h" 
#include <stdlib.h>
#include <time.h>
#include <array>
#include <fstream>

#include "boost/program_options.hpp"
namespace po = boost::program_options;





int main(int argc, char* argv[]){

 // ------------------------------ Command Line Inputs ----------------------------------------
    // Declare a group of options that will be 
    // allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
    ("version,v", "print version string")
    ("help", "produce help message")  
    ;

    // Declare a group of options that will be 
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration");
    config.add_options()
    ("Matrix_Dimension,N", po::value<int>()->default_value(5),
    "Square Matrix Dimension.")
    ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config);


    po::variables_map vm;
    store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);

    notify(vm);

    if (vm.count("help")) {

        std::cout << "Allows Options: " << std::endl;
        std::cout << generic << std::endl;
        std::cout << config << std::endl;

        return 0;
    }
    // ------------------------------- Blas Functions -------------------------------------------
    srand(time(NULL));


    const int N = vm["Matrix_Dimension"].as<int>();
    std::cout << "Chosen Dimension: "<< N << std::endl;
    double X[N*N];

    for (int i=0; i<N*N; i++) {X[i] = ((double) rand() / (RAND_MAX));}

    double norm = cblas_dnrm2(N*N, &X[0], 1);
    std::cout << "Norm: " <<norm << std::endl;
    double sum = cblas_dasum(N*N, &X[0], 1);
    std::cout << "Sum: " << sum << std::endl;

    for (int i = 0; i < N; i++){

        for (int j = 0; j < N; j++){
            std::cout << X[i*N+j] << "  ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
}