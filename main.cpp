#include <iostream>
#include <mpi/mpi.h>
#include <cassert>
#include <cmath>
#include <vector>

typedef std::pair<int, int> point;

point rank_to_coordinate(int rank)
{
	assert(rank >= 0);
	
	if (rank == 0) return point(0, 0);
	
	point x = rank_to_coordinate(rank >> 2);
	
	return point(
		(x.first  << 1) + ( rank       & 1),
		(x.second << 1) + ((rank >> 1) & 1)
	);
}

int coordinate_to_rank(point coordinate)
{
	assert(coordinate.first  >= 0);
	assert(coordinate.second >= 0);
	
	if (coordinate == point(0, 0)) return 0;
	
	int x = coordinate_to_rank(point(
		coordinate.first  >> 1,
		coordinate.second >> 1
	)) << 2;
	
	return x + (coordinate.first & 1) + ((coordinate.second & 1) << 1);
}

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
	
	
	
	// Calculate ranks of neighbors
	point coordinate = rank_to_coordinate(rank);
	
	int left = -1;
	if (coordinate.first > 0)
		left = coordinate_to_rank(point(coordinate.first - 1, coordinate.second));
	
	int right = -1;
	if (coordinate.first < np_edge - 1)
		right = coordinate_to_rank(point(coordinate.first + 1, coordinate.second));
	
	int up = -1;
	if (coordinate.second > 0)
		up = coordinate_to_rank(point(coordinate.first, coordinate.second - 1));
	
	int down = -1;
	if (coordinate.second < np_edge - 1)
		down = coordinate_to_rank(point(coordinate.first, coordinate.second + 1));
	
	std::cout
	<< "Node " << rank << " at " << rank_to_coordinate(rank)
	<< " has neighbors (" << left << ", " << right << ", " << up << ", " << down << ")"
	<< std::endl;
	
	
	
	// Grid dimensions
	assert(grid_height % np_edge == 0);
	assert(grid_width  % np_edge == 0);
	
	int local_height = grid_height / np_edge + 2;
	int local_width  = grid_width  / np_edge + 2;
	
	std::vector< std::vector<double> > local_grid(
		local_width, std::vector<double>(local_height, 0.0)
	);
	
	
	
	// Load sources
	
	
	// Loop
		// Do step
		// Communicate overlap
	
	// Output result
	
	
	MPI_Finalize();
}
