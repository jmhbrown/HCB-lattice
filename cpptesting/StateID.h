// Determines and return's "temp"'s ID.
// This allows us to index the states (columns) of our matrix.

using namespace std;

#ifndef __StateID_H

#define __StateID_H

#include <vector>
#include <stdlib.h>
#include <fstream>

// Given a vector and a matrix, returns the matrix row number which the vector matches.
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


// Makes an Nsite*Nsite identity matrix (2Darray with [width][height] = [Nsite][Nsite]
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

// Makes a row vector (an array) of lenth Nsite with zeros everywhere but loc
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


#endif
