#include <vector>

void SymmetricStored(std::vector<double> &Matr, const int n, const double alpha, const double beta);

void create_f(std::vector<double> &vec, 
                const int n, const double dx, const int a, const int b, 
                double (*func_f)(double));

std::vector<double> Conj_Grad_Method_v2(std::vector<double> &A, std::vector<double> &b, std::vector<double> &x0);


