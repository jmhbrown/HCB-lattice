#include <vector>
#include <stdlib.h>

using namespace std;

// Makes the Nsite*Nsite matrix F into an idenity matrix

void StateInit(int * F, int Nsite){
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
