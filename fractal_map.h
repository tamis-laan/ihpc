#ifndef FRACTAL_MAP_H
#define FRACTAL_MAP_H

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

#endif
