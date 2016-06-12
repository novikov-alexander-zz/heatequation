#include "mpi.h"
#include <iostream>

#include "Grid.h"
#include "Solver.h"

int rank, size;

int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int x = 10, y = 10;
	double h = 100, w = 100;
	
	x = x - 2;
	y = y - 2;

	
	Grid *myGrid = new Grid(x, y, h, w);
	Grid *tempGrid = new Grid(x, y, h, w);
	Grid *nextOne = new Grid(x, y, h, w);

	myGrid->clear();
	tempGrid->clear();
	nextOne->clear();

	Solver solver(0.10000);
	solver.solve(myGrid, nextOne, tempGrid);
	nextOne->print(0);
	MPI_Finalize();
	return 0;
}

