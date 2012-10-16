#include <vector>
#include <stdlib.h>
#include <fstream>
#include <iostream>

using namespace std;

void StateInit(int ** F, int Nsite, int diagon, int offdia) {
	for(int a=0; a<Nsite; a++) {
		for(int b=0; b<Nsite; b++) {
			if(a==b) {
				F[a][b]=diagon;
			}
			else {
				F[a][b]=offdia;
			}
		}
	}
}

// Copies Matrix B to Matrix A
void StateCopy(int ** A, int ** B, int Nsite) {
	for(int i=0; i<Nsite; i++) {
		for(int j=0; j<Nsite; j++) {
		A[i][j]=B[i][j];
		}
	}
}


int main() {
	int ** F;
	int ** N;
	int Nsite = 5;

	N=(int **)malloc(Nsite*sizeof(int *));
	for(int i = 0; i < Nsite; i++) {
	       *(N+i)=(int *)malloc(Nsite*sizeof(int));
	}

	F=(int **)malloc(Nsite*sizeof(int *));
	for(int i = 0; i < Nsite; i++) {
	       *(F+i)=(int *)malloc(Nsite*sizeof(int));
	}

	StateInit(F, Nsite,1,0);
	StateInit(N, Nsite,0,2);
	StateCopy(F, N, Nsite);

	ofstream outfile;
	outfile.open("output.txt");

	// print N
	outfile << "\tN:\n";
	for(int k = 0; k < Nsite; k++) {
		for(int l = 0; l < Nsite; l++) {
			outfile << N[k][l];
		}
		outfile << "\n";
	}

	// print F
	outfile << "\n\n\n\tF:\n";
	for(int k = 0; k < Nsite; k++) {
		for(int l = 0; l < Nsite; l++) {
			outfile << F[k][l];
		}
		outfile << "\n";
	}

	outfile.close();

	// free memory
	for(int j = 0; j < Nsite; j++) {
		free(F[j]);
	}
	free(F);

	for(int j = 0; j < Nsite; j++) {
		free(N[j]);
	}
	free(N);
	return 0;
}
