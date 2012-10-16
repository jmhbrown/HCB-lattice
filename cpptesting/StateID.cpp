// Determines and return's "temp"'s ID.
// This allows us to index the states (columns) of our matrix.

#include <vector>
#include <stdlib.h>
#include <fstream>
#include <iostream>

using namespace std;



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


int main() {
	int * temp;
	int ** N;
	int Nsite = 5;

	// Memory allocation
	temp=(int *)malloc(Nsite*sizeof(int *));
	N=(int **)malloc(Nsite*sizeof(int *));
	for(int i = 0; i < Nsite; i++) {
	       *(N+i)=(int *)malloc(Nsite*sizeof(int));
	}


	StateInit(N, Nsite);
	VectorInit(temp, Nsite, 2);

	
	ofstream outfile;
	outfile.open("output.txt");

	// print matrix N
	outfile << "\tN:\n";
	for(int k = 0; k < Nsite; k++) {
		for(int l = 0; l < Nsite; l++) {
			outfile << N[k][l];
		}
		outfile << "\n";
	}
	// print temp
	outfile<<"\n\ttemp: \n";
	for(int j = 0; j < Nsite; j++) {
		outfile << temp[j];
	}	
	// print ID
	outfile << "\n\t ID: ";
	outfile << StateId(N, temp, Nsite);



	free(temp);
	for(int j = 0; j < Nsite; j++) {
		free(N[j]);
	}
	free(N);
	return 0;

}  
