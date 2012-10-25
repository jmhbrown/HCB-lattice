// Testing how cpp deals with ranges in arrays

#include <vector>
#include <stdlib.h>
#include <lapacke.h>

using namespace std;

int main() {


	int * M;
	int * V;
	int size = 5;

	// memory allocation
	M = (int *)malloc(size*size*sizeof(int));
	V = (int *)malloc(size*sizeof(int));	


	// initialize M
	for (int a = 0; a < size; a++) {
		for (int b = 0; b < size; b++) {
			if(a!=b) {
				M[a*size+b] = 0;
			}
			else {
				M[a*size+b] = 1;
			}
		}
	}

	// initialize V
	for (int c = 0; c < size; c++) {
		V[c] = 2;
	}

	//



	// free memory
	free(M);
	free(V);

	return 0;
}
