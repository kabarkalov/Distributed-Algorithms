#include <mpi.h>
#include <iostream>
#include <vector>

#define POSITION 1

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

//Проверка топологии - можно отправлять сообщения только соседям? или нет?
int CheckTopology(MPI_Comm Comm)
{
	int ProcNum, ProcRank, RecvRank;
	MPI_Status Status;
	
	MPI_Comm_size(Comm, &ProcNum);
	MPI_Comm_rank(Comm, &ProcRank);
	if (ProcRank == 0) 
	{
		// Действия, выполняемые только процессом с рангом 0
		printf("\n Hello from process %3d", ProcRank);
		for (int i = 1; i<ProcNum; i++) 
		{
			MPI_Recv(&RecvRank, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, Comm, &Status);
			printf("\n Hello from process %3d", RecvRank);
		}
	}
	else
	{// Сообщение, отправляемое всеми процессами, кроме процесса с рангом 0
		MPI_Send(&ProcRank, 1, MPI_INT, 0, 0, Comm);
	}
	return 0;
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
		std::cout << std::endl << " ProcRank = " << ProcRank << " Neighbors["<<i<<"] = "<< Neighbors[i] <<std::endl;
	}

	return MPI_SUCCESS;

}
//Эта функция пока не написана до конца
int GetTopology(MPI_Comm Comm)
{
	int res;//Результат выполнения MPI-функций
	bool part = false;//Флаг участия в алгоритме

	int ProcId;//Номер текущего процесса
	res = MPI_Comm_rank(Comm, &ProcId);
	if (res != MPI_SUCCESS) return res;

	//Число соседей текущего процесса
	int NumNeighbors;
	res = MPI_Graph_neighbors_count(Comm, ProcId, &NumNeighbors);
	if (res != MPI_SUCCESS) return res;

	//Номера соседей текущего процесса
	int *Neighbors = new int[NumNeighbors];
	MPI_Graph_neighbors(Comm, ProcId, NumNeighbors, Neighbors);
	if (res != MPI_SUCCESS) return res;

	for (int i = 0; i < NumNeighbors; i++)
	{
		std::cout << " ProcId = " << ProcId << " Neighbors[" << i << "] = " << Neighbors[i] << std::endl;
	}

	std::vector<int> proc_known;
	proc_known.push_back(ProcId);
	std::vector<std::pair<int, int>> channels_known;

	for (int i = 0; i < NumNeighbors; i++)
	{
		channels_known.push_back({ ProcId, Neighbors[i] });
	}

	//Отправить позиции соседям

	//Обработать прием позиций от соседей
	int RecvId;
	MPI_Status Status;
	MPI_Recv(&RecvId, 1, MPI_INT, MPI_ANY_SOURCE, POSITION, Comm, &Status);

	//



	part = true;

	return MPI_SUCCESS;

}



int main(int argc, char *argv[]) 
{
	MPI_Init(&argc, &argv);
	
	MPI_Comm Comm;//Создаем свой коммуникатор со своей топологией
	int res = CreateStarTopology(MPI_COMM_WORLD,&Comm);
	if (res != MPI_SUCCESS)
	{
		std::cout <<" Something wrong while creating topology " << std::endl;
		return res;
	}
	//Проверка созданной топологии - собираем сообщения на 0-м процессе
	res = CheckTopology(Comm);
	if (res != MPI_SUCCESS)
	{
		std::cout << " Something wrong while working with topology " << std::endl;
		return res;
	}
	//Печать топологии
	res = PrintTopology(Comm);
	if (res != MPI_SUCCESS)
	{
		std::cout << " Something wrong while working with topology " << std::endl;
		return res;
	}

	MPI_Finalize();

	return 0;
}