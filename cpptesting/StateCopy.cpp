#include <vector>
#include <stdlib.h>

using namespace std;

// Copies Vector A to Vector B (both are of length Nsite)

void StateCopy(int * A, int * B, int Nsite) {
	for(int i=0; i<Nsite; i++) {
		B[i]=A[i];
	}
}
