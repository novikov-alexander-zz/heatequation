#include "mpi.h"

extern int rank, size;

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
		alpha = -alpha;///!!! чтобы соответствовать Самарскому Гулину
		double *p = new double[s_size], *q = new double[s_size];
		double *y = dst;
		double kap0 = beta/alpha, kap1 = gamma / alpha;
		double mu0 = -f[0] / alpha, mu1 = -f[(s_size - 1)*step] / alpha;
		std::cout << mu0 << " " << mu1 << std::endl;
		p[0] = kap0;
		q[0] = mu0;
		for (int i = 0; i<s_size - 2; ++i)
		{
			p[i + 1] = beta / (alpha - gamma * p[i]);
			q[i + 1] = (gamma * q[i] + f[(i+1)*step]) / (alpha - gamma * p[i]);//здесь может быть ошибка с f
		}
		y[(s_size - 1)*step] = (mu1 + kap1*q[s_size - 2]) / (1 - p[s_size - 2] * kap1);
		for (int i = s_size - 2; i >= 0; --i)
			y[i*step] = p[i] * y[(i + 1)*step] + q[i];
	}

	inline void tma(double* dst, double gamma, double alpha, double beta, int s_size, double *b){
		MPI_Status status;
		MPI_Request request;
		double tly, tx;
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
				S[i] = alpha - beta*gamma / S[i - 1];
				T[i] = (alpha*T[i - 1] - beta*gamma*t[i - 2]) / S[i - 1];
				t[i - 1] = T[i - 1] / S[i - 1];
				u[i] = T[i] - S[i] * t[i - 1];
			}
			u[1] = T[1];
			for (int i = 1; i < s_size; i++){
				u[i] = S[i] + u[i] / (u[0] + t[i - 1]);
			}
			MPI_Isend(&u[s_size - 1], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &request);
			for (int i = 0; i < s_size; i++){
				l[i] = gamma / u[i];
			}
			y[0] = b[0];
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
			for (int i = 0; i < s_size; i++){
				u[i] = S[i] + u[i] / (tu + t[i]);
			}
			for (int i = 0; i < s_size; i++){
				l[i] = gamma / u[i];
			}
			if (rank < size - 1){
				MPI_Isend(&u[s_size - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &request);
			}
			MPI_Irecv(&tly, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &request);
			y[0] = b[0] - tly;
		}

		for (int i = 0; i < s_size - 1; i++){
			y[i + 1] = b[i + 1] - l[i] * y[i];
		}

		if (rank < size - 1){
			double tly = y[s_size - 1] * l[s_size - 1];
			MPI_Send(&tly, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
			MPI_Recv(&tx, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
			x[s_size - 1] = (y[s_size - 1] - tx * beta) / u[s_size - 1];
		}
		else {
			x[s_size - 1] = y[s_size - 1] / u[s_size - 1];
		}
		for (int i = s_size - 2; i >= 0; --i){
			x[i] = (y[i] - x[i + 1] * beta) / u[i];
		}
		
		if (rank != 0)
			MPI_Isend(&x[0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &request);
		delete[] u;
		delete[] l;
		delete[] S;
		delete[] t;
		delete[] T;
		delete[] y;
		MPI_Wait(&request, &status);
	}

	inline double f(int i, int j){
		return  0; (-((double)j - 3.5)*((double)j - 3.5) + 9) / 1000;
	}

	inline double getBounds(int i, int j){
		return 10;
	}

	inline void solve(Grid *src, Grid *dst, Grid *tmp){
		double stepx = src->getXStep();
		int maxi = src->x, maxj = src->y;
		double *srcData = src->data, *dstData = dst->data, *tmpData = tmp->data;
		double hh = stepx*stepx;
		double *b = new double[maxi * maxj];
		double **bb = new double*[maxj];
		std::memset(b, 0, maxi*maxj);

		for (int j = 0; j < maxj; j++){
			bb[j] = &b[j*maxi];
		}
		for (int i = 0; i < maxi; ++i){
			for (int j = 0; j < maxj; ++j){
				b[j*maxi + i] = -srcData[j*maxi + i] / (tau / 2) - f(i, j) / 2;
			}
		}

		if (rank == 0){
			for (int j = 0; j < maxj; ++j){
				b[j*maxi] -= getBounds(0, j) / hh;
			}
		}

		if (rank == size - 1){
			for (int j = 0; j < maxj; ++j){
				b[j*maxi + maxi - 1] -= getBounds(src->x - 1, j) / hh;
			}
		}

		for (int j = 0; j < maxj; j++){
			//tma(&tmpData[j*maxi], 1.0 / hh, -1.0 / (tau / 2) - 2.0 / hh, 1.0 / hh, maxi, bb[j]);
			MPI_Barrier(MPI_COMM_WORLD);
		}

		
		std::memset(b, 0, maxi*maxj);

		for (int i = 0; i < maxi; ++i){
			for (int j = 0; j < maxj; ++j){
				b[j*maxi + i] = -tmpData[j*maxi + i] / (tau / 2) + f(i, j) / 2;
			}
		}

		for (int i = 0; i < maxi; ++i){
			b[i] -= getBounds(i, 0) / hh;
		}

		for (int i = 0; i < maxi; ++i){
			b[(maxj - 1)*maxi + i] -=  getBounds(i, src->y - 1) / hh;
		}
		


		for (int i = 0; i < maxi; i++){
			tma_seq(&dstData[i], 1.0 / hh, -1.0 / (tau / 2) - 2.0 / hh, 1.0 / hh, maxj, &b[i], maxi);
		}
		//dst->print(2);
	}

	Solver(double tau){
		this->tau = tau;
	}
};
