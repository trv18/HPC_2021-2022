#include <iostream>
using namespace std;

#include <mpi.h>

int main(int argc, char* argv[]) {


    // initialise MPI - must be called before attempting any communication
    MPI_Init(&argc, &argv);

    // Get Index of this process in the set of all process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // assigns rank to rank


    int num;
    if(rank==0){
    cout << "Enter integer " << endl;
    cin >> num;
    }

    MPI_Bcast(&num, 1, MPI_INT, 0, MPI_COMM_WORLD);

    cout << "hello from rank " << rank << " and num " << num << endl;
    // Finalise MPI process
    MPI_Finalize();

}