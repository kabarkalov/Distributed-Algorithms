#include <mpi.h>
#include <iostream>
#include <vector>

//Создание топологии "кольцо"
int CreateRingTopology(MPI_Comm InComm, MPI_Comm* RingComm)
{
	int CommSize;
	//Число процессов
	MPI_Comm_size(InComm, &CommSize);
	
	//Создаем массивы вершин и дуг для топологии "кольцо"
	int TopologySize = CommSize;//Топология охватывает все процессы
	//int index[] = { 2, 4, 6, 8, 10 };
	//Из каждой вершины выходит две дуги
	int *index = new int[TopologySize];
	for (int i = 0; i < TopologySize; i++)
	{
		index[i] = 2 * (i+1);
	}
	//	int edges[] = { 1, 4, 0, 2, 1, 3, 2, 4, 3, 0 };
	//Каждая вершина соединена с предыдущей и последующей
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
	//Создение топологии графа вида "кольцо"
	int res = MPI_Graph_create(InComm, TopologySize, index, edges, 1, RingComm);
	if (res != MPI_SUCCESS )
	{
		return res;
	}
	//Дуги и вершины больше не нужны, они запомнены в топологии
	delete[] index;
	delete[] edges;

	return MPI_SUCCESS;
}

//Создание топологии "звезда" с цетром в нулевом процессе
int CreateStarTopology(MPI_Comm InComm, MPI_Comm* StarComm)
{
	int CommSize;
	//Число процессов
	MPI_Comm_size(InComm, &CommSize);

	//Создаем массивы вершин и дуг для топологии "звезда" 
	int TopologySize = CommSize;//Топология охватывает все процессы	
	//int index[] = { 4, 5, 6, 7, 8 };
	//Первая вершина соединена со всеми остальными
	int *index = new int[TopologySize];
	index[0] = CommSize - 1;
	//Остальные вершины соединены с первой
	for (int i = 1; i < TopologySize; i++)
	{
		index[i] = CommSize - 1 + i;
	}
	//	int edges[] = {  1, 2, 3, 4, 0, 0, 0, 0  };
	//Из первой вершины выходит TopologySize - 1 дуг в остальные вершины
	int *edges = new int[(TopologySize-1)*2];
	for (int i = 0; i < TopologySize - 1; i++)
	{
		edges[i] = i + 1;
	}
	//Все остальные вершины связаны с первой
	for (int i = TopologySize - 1; i < (TopologySize - 1)*2; i++)
	{
		edges[i] = 0;
	}
	//Создение топологии графа вида "звезда"
	int res = MPI_Graph_create(InComm, TopologySize, index, edges, 1, StarComm);
	if (res != MPI_SUCCESS)
	{
		return res;
	}
	//Дуги и вершины больше не нужны, они запомнены в топологии
	delete[] index;
	delete[] edges;

	return MPI_SUCCESS;
}

//проверка топологии - можно отправлять сообщения только соседям
int CheckTopology(MPI_Comm Comm)
{

}

//Печать топологии - каждый печатает своих соседей
int PrintTopology(MPI_Comm Comm)
{
	int res;//Результат выполнения MPI-функций
	int ProcRank;//Номер текущего процесса
	res = MPI_Comm_rank(Comm, &ProcRank);	
	if (res != MPI_SUCCESS) return res;

	//Число соседей текущего процесса
	int NumNeighbors;
	res = MPI_Graph_neighbors_count(Comm, ProcRank, &NumNeighbors);
	if (res != MPI_SUCCESS) return res;

	//Номера соседей текущего процесса
	int *Neighbors = new int[NumNeighbors];
	MPI_Graph_neighbors(Comm, ProcRank,NumNeighbors, Neighbors);
	if (res != MPI_SUCCESS) return res;
	
	for (int i = 0; i < NumNeighbors; i++)
	{
		std::cout << " ProcRank = " << ProcRank << " Neighbors["<<i<<"] = "<< Neighbors[i] <<std::endl;
	}

	return MPI_SUCCESS;

}

int GetTopology(MPI_Comm Comm)
{
	int res;//Результат выполнения MPI-функций
	int part = 0;//Флаг участия в алгоритме

	int ProcRank;//Номер текущего процесса
	res = MPI_Comm_rank(Comm, &ProcRank);
	if (res != MPI_SUCCESS) return res;

	//Число соседей текущего процесса
	int NumNeighbors;
	res = MPI_Graph_neighbors_count(Comm, ProcRank, &NumNeighbors);
	if (res != MPI_SUCCESS) return res;

	//Номера соседей текущего процесса
	int *Neighbors = new int[NumNeighbors];
	MPI_Graph_neighbors(Comm, ProcRank, NumNeighbors, Neighbors);
	if (res != MPI_SUCCESS) return res;

	for (int i = 0; i < NumNeighbors; i++)
	{
		std::cout << " ProcRank = " << ProcRank << " Neighbors[" << i << "] = " << Neighbors[i] << std::endl;
	}

	std::vector<int> proc_known;
	proc_known.push_back(ProcRank);
	std::vector<std::pair<int, int>> edges_known;

	for (int i = 0; i < NumNeighbors; i++)
	{
		edges_known.push_back({ ProcRank, Neighbors[i] });
	}

	//Отправить позиции соседям

	part = 1;

	return MPI_SUCCESS;

}



int main(int argc, char *argv[]) 
{
	MPI_Init(&argc, &argv);
	
	MPI_Comm Comm;
	int res = CreateRingTopology(MPI_COMM_WORLD,&Comm);
	if (res != MPI_SUCCESS)
	{
		std::cout <<" Something wrong while creating topology " << std::endl;
		return res;
	}

	res = PrintTopology(Comm);
	if (res != MPI_SUCCESS)
	{
		std::cout << " Something wrong while printing topology " << std::endl;
		return res;
	}

	MPI_Finalize();

	return 0;
}