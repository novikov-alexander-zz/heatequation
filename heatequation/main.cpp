#include "mpi.h"
#include <iostream>

int rank, size;

class Grid {
	double step;
public:
	double *data;
	double height, width;
	int x, y;

	double getStep(){
		return step;
	}

	Grid(int x, int y, double height, double width){
		this->height = height;
		this->width = width;
		this->step = height / y;
		this->x = x;
		this->y = y;
		if (this->step != width / x){
			throw "Not uniform grid.";
		}
		data = new double[x*y];
	}
};

class Solver{
public:
	double tau;

/*	inline int solve(Grid *src, Grid *dst){
		double step = src->getStep();
		int maxi = src->x - 1, maxj = src->y - 1;
		double *srcData = src->data, *dstData = dst->data;
		int hh = step*step;
		for (int i = 1; i < maxi; ++i){
			for (int j = 1; j < maxj; ++j){
				dstData[i*maxj + j] = (/*f+*//*(srcData[i*maxj + j + 1] + srcData[i*maxj + j - 1] + srcData[(i + 1)*maxj + j] + srcData[(i - 1)*maxj + j] - 4 * srcData[i*maxj + j]) / hh)*tau + srcData[i*maxj + j];
			}
		}
		return 0;
	}
	*/
	inline void tma(double gamma, double alpha, double beta, int s_size, double *b){
		MPI_Status status;
		int portion = (s_size / (size + 1));
		int i1 = (rank == 0) ? 0 : s_size - portion * (size - rank);
		int i2 = s_size - portion * (size - 1 - rank);
		double *u = new double[i2 - i1 + 1];
		double *l = new double[i2 - i1 + 1];
		double *S = new double[i2 - i1 + 1];
		double *t = new double[i2 - i1 + 1];
		double *T = new double[i2 - i1 + 1];
		double *y = new double[i2 - i1 + 1];

		if (rank == 0){
			u[0] = alpha;
			t[0] = 0;
			S[1] = alpha;
			T[1] = -beta*gamma;
			for (int i = 2; i <= i2; i++){
				S[i] = alpha - beta*gamma/ S[i - 1];
				T[i] = (alpha*T[i - 1] - beta*gamma*t[i - 2]) / S[i - 1];
				t[i - 1] = T[i - 1] / S[i - 1];
				u[i] = T[i] - S[i] * t[i - 1];
			}
			u[1] = T[1] /*- S[1] * t[0]*/;
			for (int i = 1; i <= i2; i++){
				u[i] = S[i] + u[i] / (u[0] + t[i - 1]);
			}
			MPI_Send(&u[i2 - 1], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
			y[0] = b[0];
		}
		else {
			t[0] = 0;
			S[1] = alpha;
			T[1] = -beta*gamma;
			for (int i = 2; i <= i2 - i1; i++){
				S[i] = alpha - beta*gamma / S[i - 1];
				T[i] = (alpha*T[i - 1] - beta*gamma*t[i - 2]) / S[i - 1];
				t[i - 1] = T[i - 1] / S[i - 1];
				u[i] = T[i] - S[i] * t[i - 1];
			}
			u[1] = T[1];
			MPI_Recv(&u[0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
			std::cout << "Recvd" << u[0] << std::endl;
			for (int i = 1; i <= i2 - i1; i++){
				u[i] = S[i] + u[i] / (u[0] + t[i - 1]);
				std::cout << i + i1 << " " << u[i] << std::endl;
			}
			if (rank < size - 1){
				MPI_Send(&u[i2 - i1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
				std::cout << "Sent" << u[i2 - i1] << std::endl;
			}
		}
		for (int i = 1; i <= i2; i++){
			l[i] = gamma / u[i];
		}
	}

	inline void solve(Grid *src, Grid *dst){
		double step = src->getStep();
		int maxi = src->x - 1, maxj = src->y - 1;
		double *srcData = src->data, *dstData = dst->data;
		double hh = step*step;
		double *b = new double[maxi-2];
		tma(1.0 / hh, 1.0 / tau - 2.0 / hh, 1.0 / hh, maxi-2, b);
	}
	
	Solver(double tau){
		this->tau = tau;
	}
};

int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int x = 1000, y = 1000;
	double h = 1, w = 1;
	Solver solver(0.5);
	Grid *myGrid = new Grid(x, y, h, w);
	Grid *nextOne = new Grid(x, y, h, w);
	solver.solve(myGrid, nextOne);
	return 0;
}

