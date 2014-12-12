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

#define DEBUG 0

#define max(a,b) ((a)>(b)?a:b)

enum
{
  X_DIR, Y_DIR
};

/* global variables */
int gridsize[2];
double precision_goal;		/* precision_goal of solution */
int max_iter;			/* maximum number of iterations alowed */

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

bool write_output = true;
double omega = 1.95;		/* relaxation parameter */
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

void Debug(char *mesg, int terminate)
{
  if (DEBUG || terminate)
    printf("%d: %s\n", proc_rank, mesg);
  if (terminate)
    exit(1);
}

void Setup_Grid()
{
	int x, y, s, upper_offset[2];
	double source_x, source_y, source_val;
	FILE *f;
	
	Debug("Setup_Subgrid", 0);
	
	if (proc_rank == 0) {
		f = fopen("input.dat", "r");
		if (f == NULL)
		Debug("Error opening input.dat", 1);
		fscanf(f, "nx: %i\n", &gridsize[X_DIR]);
		fscanf(f, "ny: %i\n", &gridsize[Y_DIR]);
		fscanf(f, "precision goal: %lf\n", &precision_goal);
		fscanf(f, "max iterations: %i\n", &max_iter);
	}
	
	MPI_Bcast( gridsize,       2, MPI_INT,    0, MPI_COMM_WORLD); /* broadcast the array gridsize in one call */
	MPI_Bcast(&precision_goal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); /* broadcast precision_goal */
	MPI_Bcast(&max_iter,       1, MPI_INT,    0, MPI_COMM_WORLD); /* broadcast max_iter */
	
	/* Calculate top left corner coordinates of local grid */
	offset[X_DIR] = gridsize[X_DIR] * proc_coord[X_DIR] / P_grid[X_DIR];
	offset[Y_DIR] = gridsize[Y_DIR] * proc_coord[Y_DIR] / P_grid[Y_DIR];
	upper_offset[X_DIR] = gridsize[X_DIR] * (proc_coord[X_DIR] + 1) / P_grid[X_DIR];
	upper_offset[Y_DIR] = gridsize[Y_DIR] * (proc_coord[Y_DIR] + 1) / P_grid[Y_DIR];
	
	/* Calculate dimensions of local grid */
	dim[Y_DIR] = upper_offset[Y_DIR] - offset[Y_DIR];
	dim[X_DIR] = upper_offset[X_DIR] - offset[X_DIR];
	
	/* Add space for rows columns of neighboring grid */
	dim[Y_DIR] += 2;
	dim[X_DIR] += 2;
	
	/* allocate memory */
	if ((phi = static_cast<double **>(malloc(dim[X_DIR] * sizeof(*phi)))) == NULL)
		Debug("Setup_Subgrid : malloc(phi) failed", 1);
	if ((source = static_cast<int **>(malloc(dim[X_DIR] * sizeof(*source)))) == NULL)
		Debug("Setup_Subgrid : malloc(source) failed", 1);
	if ((phi[0] = static_cast<double *>(malloc(dim[Y_DIR] * dim[X_DIR] * sizeof(**phi)))) == NULL)
		Debug("Setup_Subgrid : malloc(*phi) failed", 1);
	if ((source[0] = static_cast<int *>(malloc(dim[Y_DIR] * dim[X_DIR] * sizeof(**source)))) == NULL)
		Debug("Setup_Subgrid : malloc(*source) failed", 1);
	
	for (x = 1; x < dim[X_DIR]; x++)
	{
		phi[x] = phi[0] + x * dim[Y_DIR];
		source[x] = source[0] + x * dim[Y_DIR];
	}
	
	/* set all values to '0' */
	for (x = 0; x < dim[X_DIR]; x++)
	for (y = 0; y < dim[Y_DIR]; y++)
	{
		phi[x][y] = 0.0;
		source[x][y] = 0;
	}
	
	/* put sources in field */
	do
	{
		if (proc_rank == 0)
			s = fscanf(f, "source: %lf %lf %lf\n", &source_x, &source_y, &source_val);
		
		MPI_Bcast(&s, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		if (s==3)
		{
			MPI_Bcast(&source_x,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&source_y,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&source_val, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			
			x = source_x * gridsize[X_DIR];
			y = source_y * gridsize[Y_DIR];
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

double Do_Step(int parity)
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
			
			delta *= omega;
			
			phi[x][y] += delta;
			
			if (max_err < fabs(delta))
				max_err = fabs(delta);
		}
	
	Exchange_Borders();
	return max_err;
}

void Solve()
{
	int count = 0;
	double delta;
	double delta1, delta2;

	Debug("Solve", 0);

	/* give global_delta a higher value then precision_goal */
	delta = 2 * precision_goal;

	while (delta > precision_goal && count < max_iter)
	{
		Debug("Do_Step 0", 0);
		delta1 = Do_Step(0);

		Debug("Do_Step 1", 0);
		delta2 = Do_Step(1);

		delta = max(delta1, delta2);
		
		MPI_Allreduce(&delta, &delta, 1, MPI_DOUBLE, MPI_MAX, grid_comm);
		
		count++;
	}

	meta << "\"iterations\" : " << count;
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
		if (argc > 3) write_output = atoi(argv[3]);
		if (argc > 4) omega        = atof(argv[4]);
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
	
	printf("(%i) (x,y)=(%i,%i)\n", proc_rank, proc_coord[X_DIR], proc_coord[Y_DIR]);
	
	/* calculate ranks of neighbouring processes */
	MPI_Cart_shift(grid_comm, X_DIR, 1, &proc_left, &proc_right);
	MPI_Cart_shift(grid_comm, Y_DIR, 1, &proc_top, &proc_bottom);
	
	if (DEBUG)
		printf("(%i) top %i, right %i, bottom %i, left %i\n",
			proc_rank, proc_top, proc_right, proc_bottom, proc_left);
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
