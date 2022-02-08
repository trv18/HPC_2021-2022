#include <iostream>
#include <vector>
#include <cmath>


#include "boost/program_options.hpp"
namespace po = boost::program_options;

#include <cblas.h>

#include "Exc_4.3_funcs.hpp"
// #include "General_funcs.hpp"

#define F77NAME(x) x##_
extern "C" {
    void F77NAME(dspev)(const char &JOBZ, const char &UPLO, // JOBZ: N  | V; UPLO: U | L
                        const int &N, double *A, double *w, // N: size; A: matrix, w: eigenvalues
                        double *Z, const int &ldz, double *WORK,  // Z: .... ; leading dim; WORK: workspace
                        int &INFO);
}

/**
 * @param [N]: amount of sample points
 * @param [dx]: sample point spacing
 * @param [U]: vector to be filled
 * @param [*func_U]: reference to function to be used
**/
void Generate_U0(const int N, const double dx, std::vector<double> &U ,double (*func_U)(double)){
    for(int i = 0; i < N; ++i){
        U[i] = func_U((i+1)*dx);
    }
}

double func_U(double x){
    return sin(M_PI*x);
}

/**
 * @brief Function to generate A matrix
 * 
 * @param N 
 * @param dx 
 * @param dt 
 * @param A : matrix of zeros to be filled in using symmetric storage. Size (N^2 + N)/2
 */
void Generate_A(const int N, const double dx, const double dt, std::vector<double> &A){
    const double nu = dt*pow(dx,-2);

    for(int i = 0; i < N ; ++i){
        for(int j = i; j < N; ++j){
            if(j==i){
                A.push_back(1-2*nu);
            }
            else if(j==i+1){
                A.push_back(nu);
            }
            else{
                A.push_back(0);
            }
        }

    }

}


int main(int argc, char* argv[]){
    typedef std::vector<double> dvec;


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
    ("Time_Dimension,T", po::value<int>()->default_value(100),
    "Time Interval to solve for.")
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

    const int N = vm["Matrix_Dimension"].as<int>();
    const int N_t = vm["Time_Dimension"].as<int>();
    bool pm = vm["Print_Matrix"].as<bool>();
    std::cout << "Chosen Dimension: "<< N << std::endl;

    // ------------------------------- Blas Functions -------------------------------------------


    dvec Matr; // coefficient matrix in Algorithm
    dvec U(N);   // matrix of U values
    dvec U_kp1(N); // matrix of U_{k+1} values
    dvec U_BC(N); // contains boundary values   
    dvec temp(N); // temporary vector to be overwritten

    const int a = 0;
    const int b = 1;

    const double dx = (b-a)/(double)(N+1); // Ensures NxN reduced matrix from (N+2)(N+2)
    const double dt = 0.001;

    Generate_U0(N, dx, U, &func_U);
    Generate_A(N, dx, dt, Matr);
    U_BC[0] = 0.0; // BC @ x = 0
    U_BC[N-1] = 0.0; // BC @ x = 1

    CBLAS_LAYOUT Layout = CblasColMajor;
    CBLAS_TRANSPOSE transA = CblasNoTrans;
    CBLAS_SIDE sideA = CblasLeft;

    const int M = N;
    const int K = N;
    const double alpha = 1.0;
    const double lda = N;
    std::vector<double> C(N*N,0) ;         // "A" matrix in Algorithm
    std::vector<double> y(N,0);    // "b" matrix in Algorithm
    const double icx = 1.0;
    const double icy = 1.0;

    temp = U_BC;
    U_kp1 = U;

    if(pm){
        print_matrix(Matr, 1, (N*N + N)/2);
        std::cout << "Original Matrix: \n";
        print_matrix(U_kp1,1,N);
    }


    for(int i=0; i<N_t; i++){
        cblas_dspmv(Layout, CblasLower, N, alpha, &Matr[0], &U_kp1[0], icx, 1.0, &temp[0], icx); // overwrites temp with U_kp1
        U_kp1 = temp;
        temp = U_BC;
    }

    if(pm){
        std::cout << "Resulting Matrix: \n";
        print_matrix(U_kp1,1,N);
    }


    // ------------------------------- LAPLAS Functions : Exc 5.2 ------------------------------------------

    dvec w(N);
    dvec Z(N);
    dvec work(3*N);
    int INFO;
    F77NAME(dspev)('N', 'L', N, &Matr[0], &w[0], &Z[0], 1, &work[0], INFO);
    std::cout << "For a value of nu = " << dt*pow(dx,-2) << " , the largest eigenvalue is: " << std::endl;    
    print_matrix(w, 1, N);
}