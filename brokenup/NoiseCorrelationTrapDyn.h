/*******************************************************************************/
#include "stdafxDyn.h"
#include <mkl_lapack.h>
#include <mkl_blas.h>
#include <fstream>
#include <sstream>
#include <cstring>


using namespace std;

/*******************Declarations of functions and global variables*****************/

void MPI_initialize(int, char **);
void StateInit(int **);
void HamGen(double *, int **, double, char);
void StateCopy(int *, int *, int);
int StateId(int **, int *);
extern void StateEvoWithTime();
extern void MomDen(cdouble *, cdouble *);
extern void NoiseCorr(cdouble *, cdouble *, cdouble *);
void output(cdouble *, cdouble *, cdouble *);
void MPI_finalize();
void EnergChk(double *, double*, double *, cdouble *, cdouble *);
void ChkOutput(double, double, double, double);
