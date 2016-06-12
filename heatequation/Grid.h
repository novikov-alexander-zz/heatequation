extern int rank, size;

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
		this->height = height;
		this->width = width;
		this->xstep = height / x;
		this->ystep = width / y;

		int portionx = (x / (size + 1));
		int i1x = (rank == 0) ? 0 : x - portionx * (size - rank);
		int i2x = x - portionx * (size - 1 - rank);
		x = i2x - i1x;

		this->x = x;
		this->y = y;
		data = new double[x*y];
	}

	void clear(){
		std::memset(data, 0, this->x*this->y*sizeof(double));
	}

	void print(){
		for (int i = 0; i < size; i++) {
			MPI_Barrier(MPI_COMM_WORLD);
			if (i == rank) {
				for (int y = 0; y < this->y; ++y){
					for (int x = 0; x < this->x; ++x){
						std::cout << data[y * this->x + x] << " ";
					}
					std::cout << std::endl << std::endl << std::endl;
				}
			}
		}
	}

	void print(int proc){
		if (proc == rank) {
			for (int y = 0; y < this->y; ++y){
				for (int x = 0; x < this->x; ++x){
					std::cout << data[y * this->x + x] << " ";
				}
				std::cout << std::endl << std::endl << std::endl;
			}
		}
	}
};