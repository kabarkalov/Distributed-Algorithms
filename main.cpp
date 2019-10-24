#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[]) 
{
	MPI_Init(&argc, &argv);
	int ProcNum, ProcRank;

	//Число процессов
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
/*
	//Создаем массивы вершин и дуг для топологии "кольцо"
	int TopologySize = ProcNum;//Топология охватывает все процессы
	//int index[] = { 2, 4, 6, 8, 10 };
	
	int *index = new int[TopologySize]; 
	for (int i = 0; i < TopologySize; i++)
	{
		index[i] = 2 * (i+1);
	}
//	int edges[] = { 1, 4, 0, 2, 1, 3, 2, 4, 3, 0 };

	int *edges = new int[2*TopologySize]; 
	edges[0] = 1;
	edges[1] = TopologySize - 1;

	for (int i = 1; i < TopologySize - 1; i++)
	{
		edges[2*i] = i - 1;
		edges[2*i+1] = i + 1;
	}
	edges[2*TopologySize - 2] = TopologySize - 2;
	edges[2*TopologySize - 1] = 0;
	
	MPI_Comm StarComm;
	int res;
	res = MPI_Graph_create(MPI_COMM_WORLD, TopologySize, index, edges, 0, &StarComm);
	if (res != MPI_SUCCESS)
	{
		std::cout <<" Something wrong while creating topology " << std::endl;
		return 1;
	}
	*/

	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	std::cout << "ProcNum = " << ProcNum << " ProcRank = " << ProcRank << std::endl;

//	delete[] index;
//	delete[] edges;

	MPI_Finalize();

	return 0;
}