#include <vector>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "MatrixUtils.h"

using namespace std;

int main() {
	double * H;
	int ** F;
	int *temp;
	int Nsite = 5;
	int i;
	
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

	return 0;
}
