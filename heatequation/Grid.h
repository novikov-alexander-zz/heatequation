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

	void Clear(){
		std::memset(data, 0, this->x*this->y*sizeof(double));
	}

	void Print(){
		for (int i = 0; i < size; i++) {
			MPI_Barrier(MPI_COMM_WORLD);
			if (i == rank) {
				for (int x = 0; x < this->x; ++x){
					for (int y = 0; y < this->y; ++y){
						std::cout << data[x * this->y + y] << " ";
					}
					std::cout << std::endl;
				} 
			}
		}
	}
	void Print(int proc){
		if (proc == rank) {
			for (int x = 0; x < this->x; ++x){
				for (int y = 0; y < this->y; ++y){
					std::cout << data[x * this->y + y];
					std::cout << " ";
				}
				std::cout << std::endl;
			}
		}
	}
		
};