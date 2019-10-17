#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[]) 
{
	int ProcNum, ProcRank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	std::cout << "ProcNum = " << ProcNum << " ProcRank = " << ProcRank << std::endl;

	MPI_Finalize();
	return 0;
}