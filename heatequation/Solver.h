#include "mpi.h"
#include <omp.h>
#include <cstring>

extern int rank, size;

class Solver{
public:
	double tau;
	double *u, *l;
	double beta;
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
	//тут f понимаетс€ в обычном смысле без минуса.
	inline void tma_seq(double *dst, double gamma, double alpha, double beta, int s_size, double *f, int step)
	{
		alpha = -alpha;///!!! чтобы соответствовать —амарскому √улину
		double *p = new double[s_size], *q = new double[s_size];
		double *y = dst;
		double kap0 = beta / alpha, kap1 = gamma / alpha;
		double mu0 = -f[0] / alpha, mu1 = -f[(s_size - 1)*step] / alpha;
		p[0] = kap0;
		q[0] = mu0;
		for (int i = 0; i<s_size - 2; ++i)
		{
			p[i + 1] = beta / (alpha - gamma * p[i]);
			q[i + 1] = (gamma * q[i] - f[(i + 1)*step]) / (alpha - gamma * p[i]);//минус перед f потому что в самарском гулине знаки наоборот!
		}
		y[(s_size - 1)*step] = (mu1 + kap1*q[s_size - 2]) / (1 - p[s_size - 2] * kap1);
		for (int i = s_size - 2; i >= 0; --i)
			y[i*step] = p[i] * y[(i + 1)*step] + q[i];
		delete[] p;
		delete[] q;
	}

	inline void tma_prepare(double gamma, double alpha, double beta, int s_size){
		MPI_Status status;
		MPI_Request request;
		if (u)
			delete[] u;
		u = new double[s_size];
		if (l)
			delete[] l;
		l = new double[s_size];
		this->beta = beta;
		double *S = new double[s_size];
		double *t = new double[s_size];
		double *T = new double[s_size];
		if (rank == 0){
			u[0] = alpha;
			t[0] = 0;
			S[1] = alpha;
			T[1] = -beta*gamma;

			for (int i = 2; i < s_size; i++){
				S[i] = alpha - beta*gamma / S[i - 1];
				T[i] = (alpha*T[i - 1] - beta*gamma*t[i - 2]) / S[i - 1];
				t[i - 1] = T[i - 1] / S[i - 1];
				u[i] = T[i] - S[i] * t[i - 1];
			}
			u[1] = T[1];
			u[s_size - 1] = S[s_size - 1] + u[s_size-1] / (u[0] + t[s_size - 2]);
			MPI_Isend(&u[s_size - 1], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &request);
#pragma omp parallel for
			for (int i = 1; i < s_size - 1; i++){
				u[i] = S[i] + u[i] / (u[0] + t[i - 1]);
				l[i] = gamma / u[i];
			}
			l[0] = gamma / u[0];
			l[s_size - 1] = gamma / u[s_size - 1];
		}
		else {
			double tu;
			MPI_Irecv(&tu, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &request);
			t[0] = 0;
			S[0] = alpha;
			T[0] = -beta*gamma;
			for (int i = 1; i < s_size; i++){
				S[i] = alpha - beta*gamma / S[i - 1];
				T[i] = (alpha*T[i - 1] - beta*gamma*t[i - 1]) / S[i - 1];
				t[i] = T[i - 1] / S[i - 1];
				u[i] = T[i] - S[i] * t[i];
			}
			u[0] = T[0];
			MPI_Wait(&request, &status);
			u[s_size - 1] = S[s_size - 1] + u[s_size - 1] / (tu + t[s_size - 1]);
			if (rank < size - 1){
				MPI_Isend(&u[s_size - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &request);
			}
#pragma omp parallel for
			for (int i = 0; i < s_size - 1; i++){
				u[i] = S[i] + u[i] / (tu + t[i]);
				l[i] = gamma / u[i];
			}
			l[s_size - 1] = gamma / u[s_size - 1];
		}
		for (int i = 1; i < s_size - 1; i++){
			l[i] = -l[i - 1] * l[i];
		}
		l[s_size - 1] = -l[s_size - 2] * l[s_size - 1];
		delete[] S;
		delete[] t;
		delete[] T;
	}

	inline void tma(double* dst, int s_size, double *b, double *um, double *y){
		MPI_Status status;
		MPI_Request srequest, rrequest;
		double tly, tx;
		int proc = omp_get_thread_num();
		double *x = dst;
		
		if (rank != 0){
			MPI_Irecv(&tly, 1, MPI_DOUBLE, rank - 1, proc, MPI_COMM_WORLD, &rrequest);
		}

		for (int i = 1; i < s_size - 1; i++){
			b[i + 1] = -l[i - 1] * b[i + 1] + b[i];
		}

		if (rank == 0){
			y[0] = b[0];
		} else {
			MPI_Wait(&rrequest, &status);
			y[0] = b[0] - tly;
		}
		/*
#pragma omp parallel{
		int threads = omp_get_num_threads();
		int thread = omp_get_thread_num();
		int portion = s_size / threads;
		int thresh = portion*(thread + 1) > s_size ? portion*(thread + 1) : s_size;
		int i = portion*thread;
		if (i == 0)
			i++;
		for (; i < thresh; ++i){
			l[i] = -l[i - 1] * l[i];
			b[i + 1] = -l[i - 1] * b[i + 1] + b[i];
		}
		int step = portion;
		portion /= 2;
		for (int i = portion*thread + 1; i < thresh; ++i){
			l[i] = -l[i - 1] * l[i];
			b[i + 1] = -l[i - 1] * b[i + 1] + b[i];
		}
	} 
*/
		y[s_size - 1] = b[s_size - 1] - l[s_size - 2] * y[0];

		if (rank < size - 1){
			double tly = y[s_size - 1] * l[s_size - 1];
			MPI_Isend(&tly, 1, MPI_DOUBLE, rank + 1, proc, MPI_COMM_WORLD, &srequest);
			MPI_Irecv(&tx, 1, MPI_DOUBLE, rank + 1, proc, MPI_COMM_WORLD, &rrequest);
		}

#pragma omp parallel for
		for (int i = 0; i < s_size - 2; i++){
			y[i + 1] = b[i + 1] - l[i] * y[0];
		}

		y[s_size - 1] = y[s_size - 1] / u[s_size - 1];
		um[s_size - 1] = -beta / u[s_size - 1];
		for (int i = s_size - 2; i >= 0; --i){
			y[i] = (-beta*y[i + 1] + y[i]) / u[i];
			um[i] = -beta*um[i + 1] / u[i];
			std::cout << um[i] << std::endl;
		}

		if (rank < size - 1){
			MPI_Wait(&rrequest, &status);
			x[s_size - 1] = y[s_size - 1] - tx * beta / u[s_size - 1];
		}
		else {
			x[s_size - 1] = y[s_size - 1];
		}

		x[0] = y[0] + um[0] * x[s_size - 1];
		if (rank != 0)
			MPI_Isend(&x[0], 1, MPI_DOUBLE, rank - 1, proc, MPI_COMM_WORLD, &srequest);

#pragma omp parallel for 
		for (int i = s_size - 2; i > 0; --i){
			x[i] = y[i] + um[i] * x[s_size - 1];
		}
	}

	inline double f(int i, int j){
		return  0; (-((double)j - 3.5)*((double)j - 3.5) + 9) / 1000;
	}

	inline double getBounds(int i, int j){
		return 10;
	}

	inline void solve(Grid *&src, Grid *&dst, Grid *tmp, int iters){
		double stepx = src->getXStep();
		int maxi = src->x, maxj = src->y;
		double *srcData = src->data, *dstData = dst->data, *tmpData = tmp->data;
		double hh = stepx*stepx;
		double *b = new double[maxi * maxj];
		double **bb = new double*[maxj];
		int threads = omp_get_max_threads();

		static double *um, *y;

#pragma omp threadprivate(um)
#pragma omp threadprivate(y)

#pragma omp parallel
		{
			um = new double[maxi];
			y = new double[maxi];
#pragma omp for
			for (int j = 0; j < maxj; j++){
				bb[j] = &b[j*maxi];
			}
		}

		

		tma_prepare(1.0 / hh, -1.0 / (tau / 2) - 2.0 / hh, 1.0 / hh, maxi);

		for (int it = 0; it < iters; ++it){

#pragma omp parallel for	
			for (int i = 0; i < maxi; ++i){
				for (int j = 0; j < maxj; ++j){
					b[j*maxi + i] = -srcData[j*maxi + i] / (tau / 2) + f(i, j) / 2;
				}
			}

			if (rank == 0){
#pragma omp parallel for	
				for (int j = 0; j < maxj; ++j){
					b[j*maxi] -= getBounds(0, j) / hh;
				}
			}

			if (rank == size - 1){
#pragma omp parallel for	
				for (int j = 0; j < maxj; ++j){
					b[j*maxi + maxi - 1] -= getBounds(src->x - 1, j) / hh;
				}
			}

#pragma omp parallel for
				for (int j = 0; j < maxj; j++){
					//tma_seq(&tmpData[j*maxi], 1.0 / hh, -1.0 / (tau / 2) 	- 2.0 / hh, 1.0 / hh, maxi, bb[j], 1);
					tma(&tmpData[j*maxi], maxi, bb[j], um, y);
				}

#pragma omp parallel for	
			for (int i = 0; i < maxi; ++i){
				for (int j = 0; j < maxj; ++j){
					b[j*maxi + i] = -tmpData[j*maxi + i] / (tau / 2) + f(i, j) / 2;
				}
				b[i] -= getBounds(i, 0) / hh;
				b[(maxj - 1)*maxi + i] -= getBounds(i, src->y - 1) / hh;

				tma_seq(&dstData[i], 1.0 / hh, -1.0 / (tau / 2) - 2.0 / hh, 1.0 / hh, maxj, &b[i], maxi);
			}
			Grid* t = dst;
			dst = src;
			src = t;
			srcData = src->data;
			dstData = dst->data;
		}
	}

	Solver(double tau){
		this->tau = tau;
		l = NULL;
		u = NULL;
	}
};
