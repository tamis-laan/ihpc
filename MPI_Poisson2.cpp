/*
 * SEQ_Poisson.c
 * 2D Poison equation solver
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi/mpi.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>
#include <algorithm>

template <typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &vec)
{
	out << "[";
	if (!vec.empty()) {
		std::copy_n(begin(vec), vec.size() - 1, std::ostream_iterator<T>(out, ","));
		out << vec.back();
	}
	return out << "]";
}

#define DEBUG 0

enum
{
  X_DIR, Y_DIR
};

/* global variables */
int gridsize[2] = {100, 100};
double precision_goal = 0.0001;	/* precision_goal of solution */
int max_iter = 5000;		/* maximum number of iterations alowed */

/* benchmark related variables */
double wtime;			/* wallclock time */

/* local grid related variables */
double **phi;			/* grid */
int **source;			/* TRUE if subgrid element is a source */
int dim[2];			/* grid dimensions */

/* process specific variables */
int proc_rank;			/* rank of current process */
int proc_coord[2];		/* coordinates of current process in processgrid */
int proc_top, proc_right, proc_bottom, proc_left; /* ranks of neigboring procs */

/* system global variables */
int P;				/* total number of processes */
int P_grid[2];			/* processgrid dimensions */
MPI_Comm grid_comm;		/* grid COMMUNICATOR */
MPI_Status status;

bool optimized_step = false;
bool write_output = true;
double omega = 1.95;		/* relaxation parameter */
int sweeps = 1;                 // number of sweeps between Exchange_Borders
std::ofstream meta;

int offset[2];
MPI_Datatype border_type[2];

void Setup_Grid();
double Do_Step(int parity);
void Solve();
void Write_Grid();
void Clean_Up();
void Debug(char *mesg, int terminate);
void start_timer();
void print_timer();

void start_timer()
{
	MPI_Barrier(MPI_COMM_WORLD);
	wtime = MPI_Wtime();
}

void print_timer()
{
	MPI_Barrier(MPI_COMM_WORLD);
	meta << ", \"time\" : " << (MPI_Wtime() - wtime);
	// ticks = clock() - ticks;
	// 100.0 * ticks * (1.0 / CLOCKS_PER_SEC) / wtime // CPU utilization
}

void Debug(const char *mesg, int terminate)
{
  if (DEBUG || terminate)
    printf("%d: %s\n", proc_rank, mesg);
  if (terminate)
    exit(1);
}

void Setup_Grid()
{
	FILE *f;
	
	Debug("Setup_Subgrid", 0);
	
	if (proc_rank == 0) {
		f = fopen("input.dat", "r");
		if (f == NULL) Debug("Error opening input.dat", 1);
	}
	
	MPI_Bcast( gridsize,       2, MPI_INT,    0, MPI_COMM_WORLD); /* broadcast the array gridsize in one call */
	MPI_Bcast(&precision_goal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); /* broadcast precision_goal */
	MPI_Bcast(&max_iter,       1, MPI_INT,    0, MPI_COMM_WORLD); /* broadcast max_iter */
	
	/* Calculate top left corner coordinates of local grid */
	offset[X_DIR] = gridsize[X_DIR] * proc_coord[X_DIR] / P_grid[X_DIR];
	offset[Y_DIR] = gridsize[Y_DIR] * proc_coord[Y_DIR] / P_grid[Y_DIR];
	
	/* Calculate dimensions of local grid */
	dim[X_DIR] = (gridsize[X_DIR] * (proc_coord[X_DIR] + 1) / P_grid[X_DIR]) - offset[X_DIR];
	dim[Y_DIR] = (gridsize[Y_DIR] * (proc_coord[Y_DIR] + 1) / P_grid[Y_DIR]) - offset[Y_DIR];
	
	/* Add space for rows columns of neighboring grid */
	dim[X_DIR] += 2;
	dim[Y_DIR] += 2;
	
	/* allocate memory */
	if ((phi = static_cast<double **>(malloc(dim[X_DIR] * sizeof(*phi)))) == NULL)
		Debug("Setup_Subgrid : malloc(phi) failed", 1);
	if ((source = static_cast<int **>(malloc(dim[X_DIR] * sizeof(*source)))) == NULL)
		Debug("Setup_Subgrid : malloc(source) failed", 1);
	if ((phi[0] = static_cast<double *>(malloc(dim[Y_DIR] * dim[X_DIR] * sizeof(**phi)))) == NULL)
		Debug("Setup_Subgrid : malloc(*phi) failed", 1);
	if ((source[0] = static_cast<int *>(malloc(dim[Y_DIR] * dim[X_DIR] * sizeof(**source)))) == NULL)
		Debug("Setup_Subgrid : malloc(*source) failed", 1);
	
	for (int x = 1; x < dim[X_DIR]; x++)
	{
		phi[x] = phi[0] + x * dim[Y_DIR];
		source[x] = source[0] + x * dim[Y_DIR];
	}
	
	/* set all values to '0' */
	for (int x = 0; x < dim[X_DIR]; x++)
	for (int y = 0; y < dim[Y_DIR]; y++)
	{
		phi[x][y] = 0.0;
		source[x][y] = 0;
	}
	
	/* put sources in field */
	int s;
	do
	{
		double source_x, source_y, source_val;
		if (proc_rank == 0)
			s = fscanf(f, "source: %lf %lf %lf\n", &source_x, &source_y, &source_val);
		
		MPI_Bcast(&s, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		if (s==3)
		{
			MPI_Bcast(&source_x,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&source_y,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&source_val, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			
			int x = source_x * gridsize[X_DIR];
			int y = source_y * gridsize[Y_DIR];
			x += 1;
			y += 1;
			
			x = x - offset[X_DIR];
			y = y - offset[Y_DIR];
			if (
				x > 0 && x < dim[X_DIR] - 1 &&
				y > 0 && y < dim[Y_DIR] - 1
			) {
				/* indices in domain of this process */
				phi[x][y] = source_val;
				source[x][y] = 1;
			}
		}
	}
	while (s==3);
	
	if (proc_rank == 0)
		fclose(f);
}

void Exchange_Borders()
{
  Debug("Exchange_Borders", 0);

  MPI_Sendrecv(&phi[1][dim[Y_DIR] - 2], 1, border_type[Y_DIR], proc_bottom, 0,
               &phi[1][             0], 1, border_type[Y_DIR], proc_top,    0,
               grid_comm, &status);

  MPI_Sendrecv(&phi[1][             1], 1, border_type[Y_DIR], proc_top,    0,
               &phi[1][dim[Y_DIR] - 1], 1, border_type[Y_DIR], proc_bottom, 0,
               grid_comm, &status);

  MPI_Sendrecv(&phi[dim[X_DIR] - 2][1], 1, border_type[X_DIR], proc_right,  0,
               &phi[             0][1], 1, border_type[X_DIR], proc_left,   0,
               grid_comm, &status);

  MPI_Sendrecv(&phi[             1][1], 1, border_type[X_DIR], proc_left,   0,
               &phi[dim[X_DIR] - 1][1], 1, border_type[X_DIR], proc_right,  0,
               grid_comm, &status);
}

double Do_Step(int parity, int n)
{
	int x, y;
	double max_err = 0.0;

	/* calculate interior of grid */
	for (x = 1; x < dim[X_DIR] - 1; x++)
	for (y = 1; y < dim[Y_DIR] - 1; y++)
		if ((offset[X_DIR] + x + offset[Y_DIR] + y) % 2 == parity && source[x][y] != 1)
		{
			double delta = (
				phi[x + 1][y] + phi[x - 1][y] +
				phi[x][y + 1] + phi[x][y - 1]
			) * 0.25 - phi[x][y];
			
			phi[x][y] += delta * omega;
			
			if (max_err < fabs(delta))
				max_err = fabs(delta);
		}
	
	if (n % sweeps == 0) Exchange_Borders();
	return max_err;
}

double Do_Step_opt(int parity, int n)
{
	double max_err = 0.0;

	int p1 = (offset[X_DIR] + offset[Y_DIR] + parity) & 1;
	for (int x = 1; x < dim[X_DIR] - 1; x++, p1 ^= 1)
	for (int y = 1 + p1; y < dim[Y_DIR] - 1; y += 2)
		if (source[x][y] != 1)
		{
			double delta = (
				phi[x + 1][y] + phi[x - 1][y] +
				phi[x][y + 1] + phi[x][y - 1]
			) * 0.25 - phi[x][y];
			
			phi[x][y] += delta * omega;
			
			if (max_err < fabs(delta))
				max_err = fabs(delta);
		}
	
	if (n % sweeps == 0) Exchange_Borders();
	return max_err;
}

void Solve()
{
	Debug("Solve", 0);

	std::vector<double> error;

	int n = 0;
	while (error.size() < max_iter)
	{
		double delta1, delta2;
		if (optimized_step) {
			delta1 = Do_Step_opt(0, ++n);
			delta2 = Do_Step_opt(1, ++n);
		} else {
			delta1 = Do_Step(0, ++n);
			delta2 = Do_Step(1, ++n);
		}
		double delta = std::max(delta1, delta2);
		MPI_Allreduce(&delta, &delta, 1, MPI_DOUBLE, MPI_MAX, grid_comm);
		error.push_back(delta);
		if (delta <= precision_goal) break;
	}

	meta << "\"iterations\" : " << error.size() << ", " << "\"errors\" : " << error;
}

void Write_Grid()
{
	std::stringstream ss;
	ss << "mpiout" << proc_rank << ".dat";
	std::ofstream file(ss.str());
	
	file << std::fixed;
	
	if (!file.good())
		Debug("Write_Grid : fopen failed", 1);
	
	Debug("Write_Grid", 0);
	
	for (int x = 1; x < dim[X_DIR] - 1; x++)
	for (int y = 1; y < dim[Y_DIR] - 1; y++)
		file
		<< (offset[X_DIR] + x) << " "
		<< (offset[Y_DIR] + y) << " "
		<< phi[x][y] << std::endl;
}

void Clean_Up()
{
	Debug("Clean_Up", 0);

	free(phi[0]);
	free(phi);
	free(source[0]);
	free(source);
}

void Setup_Proc_Grid(int argc, char **argv)
{
	int wrap_around[2];
	int reorder;
	
	Debug("My_MPI_Init", 0);
	
	/* Retrieve the number of processes P */
	MPI_Comm_size(MPI_COMM_WORLD, &P);
	
	/* Calculate the number of processes per column and per row for the grid */
	if (argc > 2)
	{
		P_grid[X_DIR] = atoi(argv[1]);
		P_grid[Y_DIR] = atoi(argv[2]);
		if (P_grid[X_DIR] * P_grid[Y_DIR] != P)
			Debug("ERROR : Proces grid dimensions do not match with P", 1);
		if (argc >  3) write_output    = atoi(argv[ 3]);
		if (argc >  4) omega           = atof(argv[ 4]);
		if (argc >  5) gridsize[X_DIR] = atoi(argv[ 5]);
		if (argc >  6) gridsize[Y_DIR] = atoi(argv[ 6]);
		if (argc >  7) precision_goal  = atof(argv[ 7]);
		if (argc >  8) max_iter        = atoi(argv[ 8]);
		if (argc >  9) sweeps          = atoi(argv[ 9]);
		if (argc > 10) optimized_step  = atoi(argv[10]);
	}
	else
		Debug("ERROR : Wrong parameterinput", 1);
	
	/* Create process topology (2D grid) */
	wrap_around[X_DIR] = 0;
	wrap_around[Y_DIR] = 0;       /* do not connect first and last process */
	reorder = 1;                  /* reorder process ranks */
	MPI_Cart_create(MPI_COMM_WORLD, 2, P_grid, wrap_around, reorder, &grid_comm);
	
	/* Retrieve new rank and Cartesian coordinates of this process */
	MPI_Comm_rank(grid_comm, &proc_rank);
	MPI_Cart_coords(grid_comm, proc_rank, 2, proc_coord);
	
	/* calculate ranks of neighbouring processes */
	MPI_Cart_shift(grid_comm, X_DIR, 1, &proc_left, &proc_right);
	MPI_Cart_shift(grid_comm, Y_DIR, 1, &proc_top, &proc_bottom);
}

void Setup_MPI_Datatypes()
{
	Debug("Setup_MPI_Datatypes", 0);
	
	/* Datatype for vertical data exchange; exchange in y-direction */
	MPI_Type_vector(dim[X_DIR] - 2, 1, dim[Y_DIR], MPI_DOUBLE, &border_type[Y_DIR]);
	MPI_Type_commit(&border_type[Y_DIR]);
	
	/* Datatype for horizontal data exchange; exchange in x-direction */
	MPI_Type_vector(dim[Y_DIR] - 2, 1, 1, MPI_DOUBLE, &border_type[X_DIR]);
	MPI_Type_commit(&border_type[X_DIR]);
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	Setup_Proc_Grid(argc, argv);
	
	if (proc_rank == 0)
		meta.open("meta.json");
	else
		meta.open("/dev/null");
	
	meta << "{";
	
	start_timer();
	Setup_Grid();
	Setup_MPI_Datatypes();
	Solve();
	if (write_output) Write_Grid();
	print_timer();
	Clean_Up();
	
	meta << "}";
	
	MPI_Finalize();
}
