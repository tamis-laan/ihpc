#include <iostream>
#include <mpi/mpi.h>
#include <cassert>
#include <cmath>
#include <vector>

template <typename T, typename U>
std::ostream &operator<<(std::ostream &out, const std::pair<T, U> &p)
{
	return out << "(" << p.first << ", " << p.second << ")";
}

int main(int argc, char **argv)
{
	// Program input
	int grid_height = 1024;
	int grid_width  = 1024;
	
	
	
	// MPI stuff
	MPI_Init(&argc, &argv);
	
	int np, rank;
	
	MPI_Comm_size(MPI_COMM_WORLD, &np);   // np should be a power of 4
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int np_edge = int(std::sqrt(np));
	
	
	
	
	
	MPI_Finalize();
}
