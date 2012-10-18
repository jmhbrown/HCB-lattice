#include <vector>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "MatrixUtils.h"

using namespace std;

int main() {
	double * H;
	int * F;
	int * temp;
	int Nsite = 5;
	int i;

	// Memory allocation
	temp=(int *)malloc(Nsite*sizeof(int));
	H=(double *)malloc(Nsite*sizeof(double));
	F=(int *)malloc(Nsite*sizeof(int));


	// initialize F and temp as identity matrix and unit vector
	StateInit(F, Nsite);
	VectorInit(temp, Nsite, 3);

	//fill H with 0s
	for(int a=0; a<Nsite; a++) {
		for(int b=0 ;b<Nsite; b++) {
			H[a*Nsite+b]=0;
		}
	}

	//Offdiagonal elements due to hopping
	for(int j=0; j<Nsite; j++) {
		for(int k = 0; k<Nsite-1; k++) {
			StateCopy(F, temp, j);
			if((F[j][k]+F[j][k+1])==1) {
				temp[k+1]=F[j][k];
				temp[k]=F[j][k+1];
				i=StateId(F, temp, Nsite);
				H[i*Nsite+j]=6.1;
			}
		}
	}

	ofstream outfile;
	outfile.open("output.txt");

	// print matrix F
	outfile << "\tF:\n";
	for(int k = 0; k < Nsite; k++) {
		for(int j = 0; j < Nsite; j++) {
			outfile << F[k][j];
		}
		outfile << "\n";
	}

	// print matrix H
	outfile << "\tH:\n";
	for(int k = 0; k < Nsite; k++) {
		for(int j = 0; j < Nsite; j++) {
			outfile << H[k*Nsite+j];
		}
		outfile << "\n";
	}

	// print temp
	outfile<<"\n\ttemp: \n";
	for(int j = 0; j < Nsite; j++) {
		outfile << temp[j];
	}



	// free memory	
	FreeIntArray(F, Nsite);
	free(H);
	free(temp);
	

	return 0;
}
