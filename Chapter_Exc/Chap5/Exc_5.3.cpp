#include <iostream>
#include <vector>
#include <cstdlib>
#include <time.h>
#include <cblas.h>

#include "boost/program_options.hpp"
namespace po = boost::program_options;

#include "boost/gil.hpp"
#include <boost/gil/extension/io/jpeg.hpp>


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
    ("n_Vecs,N", po::value<int>()->default_value(5),
    "Square Matrix Dimension.")
    ("Print_Matrix,P", po::value<bool>()->default_value(false), 
    "Output Matrices to be solved")
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

    const int k = vm["n_Vecs"].as<int>();
    bool pm = vm["Print_Matrix"].as<bool>();
    std::cout << "Using first "<< k << " vectors" << std::endl;

// ------------------------ GIL IMage processing --------------------------------

    // Read the image
    GIL::rgb8_image_t img;
    GIL::image_read_settings<GIL::jpeg_tag> read_settings;
    GIL::read_image("plane.jpg", img, read_settings);

    const int n = img.width();
    const int m = img.height();

    std::cout << "Image dimensions: " << n << " x " << m << std::endl;

    // Extract the pixel data to three matrices
    double* r = new double[n*m];
    double* g = new double[n*m];
    double* b = new double[n*m];
    auto imgv = const_view(img);
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < m; ++y) {
            r[x*m+y] = (int)GIL::at_c<0>(imgv(x,y));
            g[x*m+y] = (int)GIL::at_c<1>(imgv(x,y));
            b[x*m+y] = (int)GIL::at_c<2>(imgv(x,y));
        }
    }

    // ... do stuff here ...

    // Generate a new image to write out
    GIL::rgb8_image_t out(n, m);
    auto outv = GIL::view(out);
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < m; ++y) {
            outv(x,y) = GIL::rgb8_pixel_t(r[x*m+y], g[x*m+y], b[x*m+y]);
        }
    }
    GIL::image_write_info<GIL::jpeg_tag> write_settings;

    // GIL::write_view("output.jpg", outv, write_settings);

    delete[] r;
    delete[] g;
    delete[] b;
}

