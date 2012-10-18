// File: MatrixUtils.h
// Some utilities for playing with 2D Arrays

using namespace std;

#ifndef __MatrixUtils_H

#define __MatrixUtils_H

int StateId(int ** F, int * temp, int Nsite);
void StateInit(int ** F, int Nsite);
void VectorInit(int * V, int Nsite, int loc);
void FreeIntArray(int ** F, int Nsite);
void StateCopy(int ** F, int * temp, int row);

#endif
