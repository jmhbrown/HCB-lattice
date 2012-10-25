// A set of Utilities for working with matrices

#include <stdlib.h>
#include <vector>

using namespace std;

// Copies row of Matrix A to Matrix b (both are of length/width Nsite)
void StateCopy(int * A, int * B, int row, int Nsite) {
	for(int i=0; i<Nsite; i++) {
		for(int j=0; j<Nsite; j++) {
		A[i*Nsite+j]=B[i*Nsite+j];
		}
	}
}

// Determines and return's "temp"'s ID.
// This allows us to index the states (rows) of our matrix.
int StateId(int * F, int * temp, int Nsite) {
	bool Flag;
	int ID;
	for(ID=0; ID<Nsite; ID++) {
		Flag=true;
		for(int i=0; i<Nsite; i++) {
			if(F[ID*Nsite+i]!=temp[i]) {
				Flag=false;
				break;
			}
		}
		if(Flag) break; // Flag = true iff the IDth col of F == temp
	}
	return ID;
}

// Makes the Nsite*Nsite matrix F into an idenity matrix
void StateInit(int * F, int Nsite) {
	for(int a=0; a<Nsite; a++) {
		for(int b=0; b<Nsite; b++) {
			if(a==b) {
				F[a*Nsite+b]=1;
			}
			else {
				F[a*Nsite+b]=0;
			}
		}
	}
}
