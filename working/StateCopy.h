#include <vector>
#include <stdlib.h>

using namespace std;

// Copies Matrix B to Matrix A (both are of size Nsite*Nsite)

void StateCopy(int ** A, int ** B, int Nsite) {
	for(int i=0; i<Nsite; i++) {
		for(int j=0; j<Nsite; j++) {
		A[i][j]=B[i][j];
		}
	}
}
