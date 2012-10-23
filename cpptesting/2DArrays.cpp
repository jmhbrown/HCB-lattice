// Author : Jennifer Brown
// Version History: Made Oct 16 2012
//
//Purpose: This module is for testing out the capabilities of 2D arrays


#include <vector>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "MatrixUtils.h"

using namespace std;


int main () {
	int * temp;
	int * N;
	int Nsite = 6;

	// Memory allocation
	temp=(int *)malloc(Nsite*sizeof(int *));
	N=(int *)malloc(Nsite*Nsite*sizeof(int *));


	StateInit(N, Nsite);
	VectorInit(temp, Nsite, 2);
	
	N[1*Nsite+1]=2;
	N[2*Nsite+1]=3;
	N[3*Nsite+1]=4;

	ofstream outfile;
	outfile.open("output.txt");

	// print matrix N
	outfile << "\tN:\n";
	for(int k = 0; k < Nsite; k++) {
		for(int j = 0; j < Nsite; j++) {
			outfile << N[k*Nsite+j];
		}
		outfile << "\n";
	}


	// print temp
	outfile<<"\n\ttemp: \n";
	for(int j = 0; j < Nsite; j++) {
		outfile << temp[j];
	}

	// free memory
	free(N);	
	free(temp);
}

