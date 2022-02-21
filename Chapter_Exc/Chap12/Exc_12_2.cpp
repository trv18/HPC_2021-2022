#include <iostream>
using namespace std;

#include <mpi.h>

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    int rank;
    int numRanks;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // assigns rank to rank
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    int num;
    if(rank==0){
        cout << "Enter integer " << endl;
        cin >> num;

        cout << "Sending Message to next Rank " <<  endl;
        MPI_Ssend(&num, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
        // MPI_Recv(&num, 1, MPI_INT, numRanks-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cout << "message received at rank [" << rank+1<< "]" << endl;
    }
    else if(rank==numRanks-1){
        MPI_Recv(&num, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else{
        MPI_Recv(&num, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cout << "Sending Message to next Rank " <<  endl;
        MPI_Ssend(&num, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
        cout << "message received at rank [" << rank+1<< "]" << endl;
    }

    MPI_Finalize();

}