// File: MatrixUtils.cpp
// Some utilities for playing with 2D arrays
//

#include <vector>
#include <stdlib.h>
#include <fstream>
#include <iostream>

using namespace std;

// Copies "row" from matrix A to 1D array B. 
// Both are length Nsite
void StateCopy(int ** A, int * B, int row) {
	int Nsite = sizeof(B)/sizeof(int);
	for(int i=0; i<Nsite; i++) {
		B[i]=A[row][i];
	}
}

int StateId(int ** F, int * temp, int Nsite) {
	bool Flag;
	int ID;
	for(ID=0; ID<Nsite; ID++) {
		Flag=true;
		for(int i=0; i<Nsite; i++) {
			if(F[ID][i]!=temp[i]) {
				Flag=false;
				break;
			}
		}
		if(Flag) break; // Flag = true iff the IDth col of F == temp
	}
	return ID;
}

// makes an idenity matrix
void StateInit(int ** F, int Nsite) {
	for(int a=0; a<Nsite; a++) {
		for(int b=0; b<Nsite; b++) {
			if(a!=b) {
				F[a][b]=0;
			}
			else {
				F[a][b]=1;
			}
		}
	}
}

void VectorInit(int * temp, int Nsite, int loc) {
	for(int a = 0; a<Nsite; a++) {
		if(a!=loc) {
			temp[a] = 0;
		}
		else {
			temp[a] = 1;
		}
	}
}


void FreeIntArray(int ** M, int Nsite) {
	for(int a = 0; a<Nsite; a++) {
		free(M[a]);
	}
	free(M);
}
