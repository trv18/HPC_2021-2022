#include <vector>

void symmetric(std::vector<double> &M,  const int N);

void print_matrix(std::vector<double> &vec, const int N, const int M);

std::vector<double> Conj_Grad_Method(std::vector<double> &A, std::vector<double> &b, std::vector<double> &x0);
