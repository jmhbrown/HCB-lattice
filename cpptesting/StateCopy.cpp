#include <vector>
#include <stdlib.h>

using namespace std;

// Copies Row "row" of Matrix A to Vector B (both are of length Nsite)

void StateCopy(int ** A, int * B, int row) {
	int Nsite = sizeof(B)/sizeof(int)
	for(int i=0; i<Nsite; i++) {
		B[i]=A[row][i];
	}
}
