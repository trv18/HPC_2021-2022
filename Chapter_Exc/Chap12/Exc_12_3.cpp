#include <iostream>
#include <mpi.h>
#include <cblas.h>
#include <vector>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace std;

int main(int argc, char* argv[]) {

    // initialise MPI - must be called before attempting any communication
    MPI_Init(&argc, &argv);

    int rank;
    int numRanks;
    double dotp;
    double dotr=0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // assigns rank to rank
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    if(1024%numRanks!=0){
        cout << "Please select processes such that 1024 % n == 0 \n";
        MPI_Finalize();
        return 0;
    }
    srand(time(NULL)+rank);
    vector<double> SubArrayX(1024/numRanks, ((double) rand() / (RAND_MAX) - 0.5 )*2);
    vector<double> SubArrayY(1024/numRanks, ((double) rand() / (RAND_MAX) - 0.5 )*2);

    cout << "Rank [" << rank << "]" << " first value X: " << SubArrayX[0] << endl;
    cout << "Rank [" << rank << "]" << " first value Y: " << SubArrayY[0] << endl;

    dotp = cblas_ddot(1024/numRanks, &SubArrayX[0], 1, &SubArrayY[0], 1);
    MPI_Reduce(&dotp, &dotr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    cout << "Rank [" << rank << "]" << " value of reduced dotp: " << dotr << endl;

    MPI_Finalize();
}
