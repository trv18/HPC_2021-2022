#include <vector>

class Matrix_sum{
    private:
        int n;
        
    public:

        Matrix_sum(const int N);

        double Find_BLAS_Sum(double *matrix_value, int matrix_size);
        double Find_BLAS_Sum2(std::vector<double> &matrix);
        double Find_Manual_Sum(double (&matrix)[], int matrix_size);


};