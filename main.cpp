#include <mpi.h>
#include <iostream>
#include <vector>

//�������� ��������� "������"
int CreateRingTopology(MPI_Comm InComm, MPI_Comm* RingComm)
{
	int CommSize;
	//����� ���������
	MPI_Comm_size(InComm, &CommSize);
	
	//������� ������� ������ � ��� ��� ��������� "������"
	int TopologySize = CommSize;//��������� ���������� ��� ��������
	//int index[] = { 2, 4, 6, 8, 10 };
	//�� ������ ������� ������� ��� ����
	int *index = new int[TopologySize];
	for (int i = 0; i < TopologySize; i++)
	{
		index[i] = 2 * (i+1);
	}
	//	int edges[] = { 1, 4, 0, 2, 1, 3, 2, 4, 3, 0 };
	//������ ������� ��������� � ���������� � �����������
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
	//�������� ��������� ����� ���� "������"
	int res = MPI_Graph_create(InComm, TopologySize, index, edges, 1, RingComm);
	if (res != MPI_SUCCESS )
	{
		return res;
	}
	//���� � ������� ������ �� �����, ��� ��������� � ���������
	delete[] index;
	delete[] edges;

	return MPI_SUCCESS;
}

//�������� ��������� "������" � ������ � ������� ��������
int CreateStarTopology(MPI_Comm InComm, MPI_Comm* StarComm)
{
	int CommSize;
	//����� ���������
	MPI_Comm_size(InComm, &CommSize);

	//������� ������� ������ � ��� ��� ��������� "������" 
	int TopologySize = CommSize;//��������� ���������� ��� ��������	
	//int index[] = { 4, 5, 6, 7, 8 };
	//������ ������� ��������� �� ����� ����������
	int *index = new int[TopologySize];
	index[0] = CommSize - 1;
	//��������� ������� ��������� � ������
	for (int i = 1; i < TopologySize; i++)
	{
		index[i] = CommSize - 1 + i;
	}
	//	int edges[] = {  1, 2, 3, 4, 0, 0, 0, 0  };
	//�� ������ ������� ������� TopologySize - 1 ��� � ��������� �������
	int *edges = new int[(TopologySize-1)*2];
	for (int i = 0; i < TopologySize - 1; i++)
	{
		edges[i] = i + 1;
	}
	//��� ��������� ������� ������� � ������
	for (int i = TopologySize - 1; i < (TopologySize - 1)*2; i++)
	{
		edges[i] = 0;
	}
	//�������� ��������� ����� ���� "������"
	int res = MPI_Graph_create(InComm, TopologySize, index, edges, 1, StarComm);
	if (res != MPI_SUCCESS)
	{
		return res;
	}
	//���� � ������� ������ �� �����, ��� ��������� � ���������
	delete[] index;
	delete[] edges;

	return MPI_SUCCESS;
}

//�������� ��������� - ����� ���������� ��������� ������ �������
int CheckTopology(MPI_Comm Comm)
{

}

//������ ��������� - ������ �������� ����� �������
int PrintTopology(MPI_Comm Comm)
{
	int res;//��������� ���������� MPI-�������
	int ProcRank;//����� �������� ��������
	res = MPI_Comm_rank(Comm, &ProcRank);	
	if (res != MPI_SUCCESS) return res;

	//����� ������� �������� ��������
	int NumNeighbors;
	res = MPI_Graph_neighbors_count(Comm, ProcRank, &NumNeighbors);
	if (res != MPI_SUCCESS) return res;

	//������ ������� �������� ��������
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
	int res;//��������� ���������� MPI-�������
	int part = 0;//���� ������� � ���������

	int ProcRank;//����� �������� ��������
	res = MPI_Comm_rank(Comm, &ProcRank);
	if (res != MPI_SUCCESS) return res;

	//����� ������� �������� ��������
	int NumNeighbors;
	res = MPI_Graph_neighbors_count(Comm, ProcRank, &NumNeighbors);
	if (res != MPI_SUCCESS) return res;

	//������ ������� �������� ��������
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

	//��������� ������� �������

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