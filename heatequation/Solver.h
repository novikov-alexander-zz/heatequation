#include "mpi.h"
#include <omp.h>
#include <cstring>

extern int rank, size;

class Solver{
public:
	double tau;
	double *u, *l;
	double leftl;
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
	//тут f понимается в обычном смысле без минуса.
	inline void tma_seq(double *dst, double gamma, double alpha, double beta, int s_size, double *f, int step, double* p, double* q)
	{
		alpha = -alpha;///!!! чтобы соответствовать Самарскому Гулину
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
	}

	inline void tma_prepare(double gamma, double alpha, double beta, int tx, int s_size){
		MPI_Status status;
		MPI_Request rrequest, srequest;
		if (u)
			delete[] u;
		u = new double[s_size];
		if (l)
			delete[] l;
		l = new double[s_size];
		this->beta = beta;
#pragma omp parallel for
			for (int i = 0; i < s_size; i++){
				u[i] = alpha;
			}
			leftl = pow(-gamma / alpha, tx);
			l[0] = -pow( - gamma / alpha, tx+1);
		for (int i = 1; i < s_size; i++){
			l[i] = -l[i - 1] * l[i];//тут все верно
		}
	}

	inline void tma(double* dst, int s_size, double *b, double *um, double *y){
		MPI_Status status;
		MPI_Request srequest, rrequest;
		double y_to_r, y_from_l, b_from_l, b_to_r, tx;
		//double tly, tx,stly, rmul, ll, lb;
		int proc = omp_get_thread_num();
		double *x = dst;

		int owners = 1;

		if (rank == 0){//делаем его новичком
			for (int i = 0; i < s_size - 1; i++){
				b[i + 1] += -l[i] * b[i] ;
			}

			y[0] = b[0];

#pragma omp parallel for
			for (int i = 0; i < s_size - 1; i++){
				y[i + 1] = b[i + 1] - l[i] * y[0];
			} //считает себе

			y_to_r = y[s_size - 1];
		}
		else {

			if (rank == 1){//готовим принимать y
				for (int i = 1; i < s_size - 1; i++){
					b[i + 1] += -l[i] * b[i];
				}
			}

			while (rank > owners){
				if (rank % (2 * owners) < owners){//эти отправляют b
					if (rank / owners < owners / 2){//старички
						MPI_Recv(&b_to_r, 1, MPI_DOUBLE, rank + owners / 2, proc, MPI_COMM_WORLD, &status);
						MPI_Send(&b_to_r, 1, MPI_DOUBLE, rank + owners, proc, MPI_COMM_WORLD);
					}
					else {//новички
						MPI_Send(&b[s_size - 1], 1, MPI_DOUBLE, rank + owners, proc, MPI_COMM_WORLD);
					}
				}
				else {//эти принимают
					if (rank < owners * 2){//принимают y
						MPI_Recv(&y_from_l, 1, MPI_DOUBLE, rank - owners, proc, MPI_COMM_WORLD, &status);

						y[0] = b[0] - leftl*y_from_l;
						y[s_size - 1] = b[s_size - 1] - l[s_size - 2] * y_from_l;

						if (owners * 2 < size){//посылает, что получилось приславшему
							y_to_r = y[s_size - 1];
							MPI_Isend(&y_to_r, 1, MPI_DOUBLE, rank - owners, proc, MPI_COMM_WORLD, &srequest);
						}

#pragma omp parallel for
						for (int i = 0; i < s_size - 2; i++){
							y[i + 1] = b[i + 1] - l[i] * y[0];
						} //считает себе

						if (owners * 2 < size){
							MPI_Wait(&srequest, &status);
						}
					}
					else {//принимают b
						MPI_Recv(&b_from_l, 1, MPI_DOUBLE, rank - owners, proc, MPI_COMM_WORLD, &status);
						b[0] += -leftl*b_from_l;
						for (int i = 1; i < s_size - 1; i++){
							b[i + 1] += -l[i] * b[i];
						}
						if (owners * 2 < size){//посылает, что получилось приславшему
							MPI_Send(&b[s_size - 1], 1, MPI_DOUBLE, rank - owners, proc, MPI_COMM_WORLD);
						}
					}
					owners *= 2;
				}
			}
		}
		
		if (rank < owners){//шлет свой y
			MPI_Send(&y_to_r, 1, MPI_DOUBLE, rank + owners, proc, MPI_COMM_WORLD);
			owners *= 2;
		}

		while (owners < size){
			MPI_Recv(&y_to_r, 1, MPI_DOUBLE, rank + owners / 2, proc, MPI_COMM_WORLD, &status);
			MPI_Send(&y_to_r, 1, MPI_DOUBLE, rank + owners, proc, MPI_COMM_WORLD);
			owners *= 2;
		}
		


		//переходим к обратному шагу

		if (rank == size - 1){//делаем его новичком
			y[s_size - 1] = y[s_size - 1] / u[s_size - 1];
			um[s_size - 1] = -beta / u[s_size - 1];
			for (int i = s_size - 2; i >= 0; --i){
				y[i] = (-beta*y[i + 1] + y[i]) / u[i];
				um[i] = -beta*um[i + 1] / u[i];
			}
			x[s_size - 1] = y[s_size - 1];
			x[0] = y[0] + um[0] * x[s_size - 1];

#pragma omp parallel for 
			for (int i = s_size - 2; i > 0; --i){
				x[i] = y[i] + um[i] * x[s_size - 1];
			}
			//считает себе
		}
		if (rank == 1){
			y[s_size - 1] = y[s_size - 1] / u[s_size - 1];
			um[s_size - 1] = -beta / u[s_size - 1];
			for (int i = s_size - 2; i >= 0; --i){
				y[i] = (-beta*y[i + 1] + y[i]) / u[i];
				um[i] = -beta*um[i + 1] / u[i];
			}
		}

		int sizem = size - 1;
		for (int owners = 1; owners < size; owners *= 2){
			if (sizem - rank < owners){//они шлют x
				if (sizem - rank < owners / 2){//старички
					MPI_Recv(&x[0], 1, MPI_DOUBLE, rank - owners / 2, proc, MPI_COMM_WORLD, &status);
					MPI_Send(&x[0], 1, MPI_DOUBLE, rank - owners, proc, MPI_COMM_WORLD);
				}
				else {//новички
					MPI_Send(&x[0], 1, MPI_DOUBLE, rank - owners, proc, MPI_COMM_WORLD);
				}
			}
			else if ((sizem - rank) % (2 * owners) < owners){//эти отправляют um
				if (sizem - rank < owners / 2){//старички
					MPI_Recv(&um[0], 1, MPI_DOUBLE, rank - owners / 2, proc, MPI_COMM_WORLD, &status);
					MPI_Send(&um[0], 1, MPI_DOUBLE, rank - owners, proc, MPI_COMM_WORLD);
				}
				else {//новички
					MPI_Send(&um[0], 1, MPI_DOUBLE, rank - owners, proc, MPI_COMM_WORLD);
				}
			}
			else {//эти принимают
				if (sizem - rank < owners * 2){//принимают x
					MPI_Recv(&tx, 1, MPI_DOUBLE, rank + owners, proc, MPI_COMM_WORLD, &status);

					x[s_size - 1] = y[s_size - 1] - tx * beta / u[s_size - 1];
					x[0] = y[0] + um[0] * x[s_size - 1];

					if (owners * 2 < size){
						MPI_Isend(&x[0], 1, MPI_DOUBLE, rank + owners, proc, MPI_COMM_WORLD, &srequest);
					}
#pragma omp parallel for 
					for (int i = s_size - 2; i > 0; --i){
						x[i] = y[i] + um[i] * x[s_size - 1];
					}//считает себе
					if (owners * 2 < size){
						MPI_Wait(&srequest, &status);
					}
				}
				else {//принимают um
					MPI_Recv(&um[s_size - 1], 1, MPI_DOUBLE, rank + owners, proc, MPI_COMM_WORLD, &status);
					y[s_size - 1] = y[s_size - 1] / u[s_size - 1];
					um[s_size - 1] = -beta / u[s_size - 1];
					for (int i = s_size - 2; i >= 0; --i){
						y[i] = (-beta*y[i + 1] + y[i]) / u[i];
						um[i] = -beta*um[i + 1] / u[i];
					}
					if (owners * 2 < size){//посылает, что получилось приславшему
						MPI_Send(&um[s_size - 1], 1, MPI_DOUBLE, rank - owners, proc, MPI_COMM_WORLD);
					}
				}
			}
		}
	}

	inline double f(int i, int j){
		return  0; (-((double)j - 3.5)*((double)j - 3.5) + 9) / 1000;
	}

	inline double getBounds(int i, int j){
		return 10;
	}

	inline void solve(Grid *&src, Grid *&dst, Grid *tmp, int iters){
		double stepx = src->getXStep(), stepy = src->getYStep();
		int maxi = src->x, maxj = src->y;
		double *srcData = src->data, *dstData = dst->data, *tmpData = tmp->data;
		double xx = stepx*stepx, yy = stepy*stepy;
		double *b = new double[maxi * maxj];
		double **bb = new double*[maxj];
		int threads = omp_get_max_threads();

		static double *um, *y, *p, *q;

#pragma omp threadprivate(um)
#pragma omp threadprivate(y)
#pragma omp threadprivate(p)
#pragma omp threadprivate(q)

#pragma omp parallel
		{
			um = new double[maxi];
			y = new double[maxi];
			p = new double[maxj];
			q = new double[maxj];
#pragma omp for
			for (int j = 0; j < maxj; j++){
				bb[j] = &b[j*maxi];
			}
		}

		tma_prepare(1.0 / xx, -1.0 / (tau / 2) - 2.0 / xx, 1.0 / xx, src->overallx, maxi);

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
					b[j*maxi] -= getBounds(0, j) / xx;
				}
			}

			if (rank == size - 1){
#pragma omp parallel for	
				for (int j = 0; j < maxj; ++j){
					b[j*maxi + maxi - 1] -= getBounds(src->x - 1, j) / xx;
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
				b[i] -= getBounds(i, 0) / yy;
				b[(maxj - 1)*maxi + i] -= getBounds(i, src->y - 1) / yy;

				tma_seq(&dstData[i], 1.0 / yy, -1.0 / (tau / 2) - 2.0 / yy, 1.0 / yy, maxj, &b[i], maxi, p, q);
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
