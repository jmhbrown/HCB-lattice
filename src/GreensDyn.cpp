/*	The two point green's function and the four point green's function for spins are implemented here	*/

#include "stdafxDyn.h"
#include <mkl_lapack.h>
#include <mkl_blas.h>

using namespace std;


cdouble TwoPntGreen(int, int);

/*******************Function definitions***************************/
void StateEvoWithTime() //calculate the time-dependent fermionic slater determinant
{
	//P_t=U*e^(-iDt/h)*(U+)*P_0; slater determinant in time evolution
	double * P0, * Ptemp;
	cdouble * cHamT, *cPtemp;

	P0=(double *)malloc(Nfermion*Nsite*sizeof(double));
	Ptemp=(double *)malloc(Nfermion*Nsite*sizeof(double));

	cHamT=(cdouble *)malloc(Nsite*Nsite*sizeof(cdouble));
	cPtemp=(cdouble *)malloc(Nfermion*Nsite*sizeof(cdouble));
	cPT=(cdouble *)malloc(Nfermion*Nsite*sizeof(cdouble));

	//P_0 the initial fermionic slater determinant
	for(int i=0; i<Nfermion*Nsite; i++)
		P0[i]=Ham0[i];

	//Ptemp=(U+)*P_0
	double alpha=1.0, beta=0.0;
	dgemm("C", "N", &Nsite, &Nfermion, &Nsite, &alpha, HamT, &Nsite, P0, &Nsite, &beta, Ptemp, &Nsite);

	//cPtemp=1e^(-iDt)*Ptemp
	for(int i=0; i<Nfermion; i++)
		for(int j=0; j<Nsite; j++)
			cPtemp[i*Nsite+j]=polar(1.0, -double(xtime)*EnergT[j])*Ptemp[i*Nsite+j];

	//P_t=U*cPtemp
	for(int i=0; i<Nsite*Nsite; i++)
		cHamT[i]=HamT[i];
	cdouble calpha=1.0, cbeta=0.0;
	zgemm("N", "N", &Nsite, &Nfermion, &Nsite, &calpha, cHamT, &Nsite, cPtemp, &Nsite, &cbeta, cPT, &Nsite);

	free(P0);
	free(Ptemp);
	free(cHamT);
	free(cPtemp);
}

cdouble TwoPntGreen(int i, int j)
{
	cdouble * Pa, * Pb; //They represents states revised from P, with dimension (Nfermion+1)*Nsite
	cdouble * C; //Store the product of Pa'*Pb
	cdouble rho=1.0, alpha, beta;
	int m, n, k, lda, ldb, ldc, info, *piv; //<------for Lapack use

	Pa=(cdouble *)malloc((Nfermion+1)*Nsite*sizeof(cdouble));
	Pb=(cdouble *)malloc((Nfermion+1)*Nsite*sizeof(cdouble));
	C=(cdouble *)malloc((Nfermion+1)*(Nfermion+1)*sizeof(cdouble));
	piv=(int *)malloc((Nfermion+1)*sizeof(int));
	
	//Construct Pa and Pb by revising the multi-fermion counterpart P
	//The revision includes sign-changing for part of the elements and an addition of one column
	for(int k1=0; k1<Nfermion; k1++) 
	{
		int k2=0;
		while(k2<i)
		{
			Pa[k1*Nsite+k2]=-cPT[k1*Nsite+k2]; //sign changes
			k2++;
		}
		for(;k2<Nsite; k2++)
			Pa[k1*Nsite+k2]=cPT[k1*Nsite+k2];
	}

	for(int k1=0; k1<Nfermion; k1++)
	{
		int k2=0;
		while(k2<j)
		{
			Pb[k1*Nsite+k2]=-cPT[k1*Nsite+k2]; //sign changes
			k2++;
		}
		for(;k2<Nsite; k2++)
			Pb[k1*Nsite+k2]=cPT[k1*Nsite+k2];
	}

	for(int k2=0; k2<Nsite; k2++)
	{
		if(k2==i)
			Pa[Nsite*Nfermion+k2]=1.0;
		else
			Pa[Nsite*Nfermion+k2]=0.0;
		if(k2==j)
			Pb[Nsite*Nfermion+k2]=1.0;
		else
			Pb[Nsite*Nfermion+k2]=0.0;
	}

	//Calling lapack's cgemm() to do the matrix product
	ldc=m=n=Nfermion+1;
	lda=ldb=k=Nsite;
	alpha=1.0;
	beta=0.0;

	zgemm("C", "N", &m, &n, &k, &alpha, Pa, &lda, Pb, &ldb, &beta, C, &ldc);

	zgetrf(&m, &n, C, &ldc, piv, &info);

	if(info!=0)
		rho=0.0;
	else
		for(int h=0; h<(Nfermion+1);h++)
		{
			if((piv[h]-1)!=h)	
				rho=-rho*C[h*(Nfermion+1)+h];
			else
				rho=rho*C[h*(Nfermion+1)+h];
		}

	free(Pa);
	free(Pb);
	free(C);
	free(piv);

	return rho;
}


void TwoPntGreenMx(cdouble * g_2) //green_2 is Nsite*Nsite in dimension
{
	int flag=0;

	for(int i=0; i<Nsite; i++)
		for(int j=0; j<Nsite; j++)
			g_2[i*Nsite+j]=0.0;

	for(int i=0; i<Nsite; i++)
		for(int j=i; j<Nsite; j++)
		{
			if(flag%isize==irank)		//<------------ when the job is parallelized, this is part of MPI to recongize each CPU.
			{
				g_2[i*Nsite+j]=TwoPntGreen(i, j);
				if(i!=j)
					g_2[j*Nsite+i]=conj(g_2[i*Nsite+j]);
			}

		#ifdef MYMPI

			if(irank==0&&flag%isize!=irank)
			{
                MPI::COMM_WORLD.Recv(&g_2[i*Nsite+j], 1, MPI::DOUBLE_COMPLEX, flag%isize, i*Nsite+j, istatus);
				if(j!=i)
					g_2[j*Nsite+i]=conj(g_2[i*Nsite+j]);
			}
			else if(irank!=0&&flag%isize==irank)
				MPI::COMM_WORLD.Send(&g_2[i*Nsite+j], 1, MPI::DOUBLE_COMPLEX, 0, i*Nsite+j);
		#endif

			flag++;
		}

#ifdef MYMPI
	MPI::COMM_WORLD.Bcast(g_2, Nsite*Nsite, MPI::DOUBLE_COMPLEX, 0);
#endif
//	cout<<"rank"<<irank<<": the flag of rho matrix is "<<flag<<endl;
}


/************************************************************En*d*******************************************************************/
