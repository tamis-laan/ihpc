#include <iostream>
#include <mpi/mpi.h>

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	
	int np, rank;
	
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	std::cout << "Node " << rank << " of " << np << " says: Hello world!" << std::endl;
	
	MPI_Finalize();
	return 0;
}
