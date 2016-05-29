#include "mpi.h"
#include <iostream>

int rank, size;

class Grid {
	double xstep, ystep;
public:
	double *data;
	double height, width;
	int x, y;

	double getXStep(){
		return xstep;
	}

	double getYStep(){
		return xstep;
	}

	Grid(int x, int y, double height, double width){
		int portionx = (x / (size + 1));
		int i1x = (rank == 0) ? 0 : x - portionx * (size - rank);
		int i2x = x - portionx * (size - 1 - rank);
		x = i2x - i1x;

		this->height = height;
		this->width = width;
		this->xstep = height / x;
		this->ystep = height / y;
		this->x = x;
		this->y = y;
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

	inline void tma_seq(double *dst, double gamma, double alpha, double beta, int s_size, double *f, int step)
	{
		double *p = new double[s_size], *q = new double[s_size];
		double *y = dst;
		p[0] = 0;
		q[0] = f[0];
		for (int i = 0; i<s_size-2; ++i)
		{
			p[i + 1] = beta / (alpha - gamma * p[i]);
			q[i + 1] = (gamma * q[i] + f[i*step]) / (alpha - gamma * p[i]);
		}
		y[(s_size - 1)*step] = f[(s_size - 1)*step];
		for (int i = s_size - 2; i >= 0; --i)
			y[i*step] = p[i] * y[(i + 1)*step] + q[i];
	}

	inline void tma(double* dst, double gamma, double alpha, double beta, int s_size, double *bá){
		MPI_Status status;
		double ty,tx;
		double *u = new double[s_size];
		double *l = new double[s_size];
		double *S = new double[s_size];
		double *t = new double[s_size];
		double *T = new double[s_size];
		double *y = new double[s_size];
		double *x = dst;
		
		if (rank == 0){
			u[0] = alpha;
			t[0] = 0;
			S[1] = alpha;
			T[1] = -beta*gamma;

			for (int i = 2; i < s_size; i++){
				S[i] = alpha - beta*gamma/ S[i - 1];
				T[i] = (alpha*T[i - 1] - beta*gamma*t[i - 2]) / S[i - 1];
				t[i - 1] = T[i - 1] / S[i - 1];
				u[i] = T[i] - S[i] * t[i - 1];
			}
			u[1] = T[1];
			for (int i = 1; i < s_size; i++){
				u[i] = S[i] + u[i] / (u[0] + t[i - 1]);
			}
			MPI_Send(&u[s_size - 1], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
			for (int i = 1; i < s_size; i++){
				l[i] = gamma / u[i];
			}
			y[0] = b[0];
		}
		else {
			double tu;
			t[0] = 0;
			S[0] = alpha;
			T[0] = -beta*gamma;
			for (int i = 1; i < s_size; i++){
				S[i] = alpha - beta*gamma / S[i - 1];
				T[i] = (alpha*T[i - 1] - beta*gamma*t[i - 2]) / S[i - 1];
				t[i] = T[i - 1] / S[i - 1];
				u[i] = T[i] - S[i] * t[i];
			}
			u[0] = T[0];
			
			MPI_Recv(&tu, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
			for (int i = 0; i < s_size; i++){
				u[i] = S[i] + u[i] / (tu + t[i]);

			}
			if (rank < size - 1){
				MPI_Send(&u[s_size - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
			}
			for (int i = 0; i < s_size; i++){
				l[i] = gamma / u[i];
			}
			MPI_Recv(&ty, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
			y[0] = b[0] - l[0] * ty;
		}

		for (int i = 0; i < s_size - 1; i++){
			y[i + 1] = b[i + 1] - l[i + 1] * y[i];
		}
		
		if (rank < size - 1){
			MPI_Send(&y[s_size - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
			MPI_Recv(&tx, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
			x[s_size - 1] = (y[s_size - 1] - tx * beta) / u[s_size - 1];
		} else {
			x[s_size - 1] = y[s_size - 1] / u[s_size - 1];
		}
		for (int i = s_size - 2; i >= 0; --i){
			x[i] = (y[i] - x[i + 1] * beta) / u[i];
		}
		if (rank!=0)
			MPI_Send(&x[0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
		delete[] u;
		delete[] l;
		delete[] S;
		delete[] t;
		delete[] T;
		delete[] y;
	}

	inline double f(int i, int j){
		return 0;
	}

	inline double getBounds(int i, int j){
		return 10;
	}

	inline void solve(Grid *src, Grid *dst, Grid *tmp){
		double stepx = src->getXStep();
		int maxi = src->x - 1, maxj = src->y - 1;
		double *srcData = src->data, *dstData = dst->data, *tmpData = tmp->data;
		double hh = stepx*stepx;
		double *b = new double[maxi * maxj];
		double **bb = new double*[maxj];
		for (int j = 0; j < maxi*maxj; j++){
			b[j] = 0;
		}
		for (int j = 0; j < maxj; j++){
			bb[j] = &b[j*maxi];
		}
		for (int i = 0; i < maxi; ++i){
			for (int j = 0; j < maxj; ++j){
				b[i] = -srcData[i*maxj + j] / (tau / 2) - f(i, j) / 2;
			}
		}
		
		if (rank == 0){
			for (int i = 0; i < maxi; ++i){
				b[i] -= getBounds(0, 1 + i) / hh;
			}
		}

		if (rank == size - 1){
			for (int i = 0; i < maxi; ++i){
				b[i + maxj - 1] -= getBounds(src->x - 1, 1 + i) / hh;
			}
		}

		for (int j = 0; j < maxj; j++){
			tma(tmpData, 1.0 / hh, 1.0 / (tau/2) - 2.0 / hh, 1.0 / hh, maxi - 1, bb[j]);
			MPI_Barrier(MPI_COMM_WORLD);
		}

		for (int i = 0; i < maxi; i++){
			tma_seq(1.0 / hh, 1.0 / (tau/2) - 2.0 / hh, 1.0 / hh, maxj, &b[i], maxi);
		}
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
	
	x = x - 2;
	y = y - 2;

	
	Grid *myGrid = new Grid(x, y, h, w);
	Grid *tempGrid = new Grid(x, y, h, w);
	Grid *nextOne = new Grid(x, y, h, w);

	Solver solver(0.5);
	solver.solve(myGrid, nextOne, tempGrid);
	MPI_Finalize();
	return 0;
}

