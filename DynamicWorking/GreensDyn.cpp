/*	The two point green's function and the four point green's 
 *	function for spins are implemented here	*/

#include "stdafxDyn.h"
#include <mkl_lapack.h>
#include <mkl_blas.h>

using namespace std;


cdouble TwoPntGreen(int, int);
cdouble FourPntGreen(int, int, int, int);

/***Function definitions***/



cdouble TwoPntGreen(int i, int j) // calculating the green's function from eq. (6.12) in Rigol 2004 -- J 
{
	cdouble * Pa, * Pb; //They represents states revised from P, with dimension (Nfermion+1)*Nsite
	cdouble * C; //Store the product of Pa'*Pb
	cdouble rho=1.0, alpha, beta;
	// Lapack use
	int m, n, k, lda, ldb, ldc, info;
	int *piv; // keeps track of row changes for our decomposition

	Pa=(cdouble *)malloc((Nfermion+1)*Nsite*sizeof(cdouble));
	Pb=(cdouble *)malloc((Nfermion+1)*Nsite*sizeof(cdouble));
	C=(cdouble *)malloc((Nfermion+1)*(Nfermion+1)*sizeof(cdouble));
	piv=(int *)malloc((Nfermion+1)*sizeof(int));
	
	//Construct Pa and Pb by revising the multi-fermion counterpart P
	//The revision includes sign-changing for part of the elements and an addition of one column
	//
	//This is converting the multiparticle state equation for fermions
	// into something applicable to HCB. See (6.10) and (6.11) --J
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

// Adding an extra column, which is motivated by the fi+ near the end of (6.6)
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

	//Calling lapack's zgemm() to do the matrix product
	ldc=m=n=Nfermion+1;
	lda=ldb=k=Nsite;
	alpha=1.0;
	beta=0.0;

	zgemm("C", "N", &m, &n, &k, &alpha, Pa, &lda, Pb, &ldb, &beta, C, &ldc);

	zgetrf(&m, &n, C, &ldc, piv, &info); // does the domposition, so we can get the determinant

	if(info!=0)   // info gets information about errors - J
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
} // two point green's function


// Not useful to Jenny. 

cdouble FourPntGreen(int i, int j,int a, int b)
{
	cdouble * Pa, * Pb; //They represents states revised from P, with dimension (Nfermion+2)*Nsite
	cdouble * C; //Store the product of Pa'*Pb
	cdouble rho=1.0, alpha, beta;
	int m, n, k, lda, ldb, ldc, info, *piv;

	Pa=(cdouble *)malloc((Nfermion+2)*Nsite*sizeof(cdouble));
	Pb=(cdouble *)malloc((Nfermion+2)*Nsite*sizeof(cdouble));
	C=(cdouble *)malloc((Nfermion+2)*(Nfermion+2)*sizeof(cdouble));
	piv=(int *)malloc((Nfermion+2)*sizeof(int));
	
	//Construct Pa and Pb by revising the multi-fermion counterpart P
	//The revision includes sign-changing for part of the elements and an addition of two columns
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
		while(k2<b)
		{
			Pb[k1*Nsite+k2]=-cPT[k1*Nsite+k2]; //sign changes
			k2++;
		}
		for(;k2<Nsite; k2++)
			Pb[k1*Nsite+k2]=cPT[k1*Nsite+k2];
	}

	for(int k2=0; k2<Nsite; k2++)	//add the first column
	{
		if(k2==i)
			Pa[Nsite*Nfermion+k2]=1.0;
		else
			Pa[Nsite*Nfermion+k2]=0.0;
		if(k2==b)
			Pb[Nsite*Nfermion+k2]=1.0;
		else
			Pb[Nsite*Nfermion+k2]=0.0;
	}
	
	//second step of revision
	for(int k1=0; k1<Nfermion+1; k1++) 
	{
		int k2=0;
		while(k2<j)
		{
			Pa[k1*Nsite+k2]=-Pa[k1*Nsite+k2]; //sign changes
			k2++;
		}
	}

	for(int k1=0; k1<Nfermion+1; k1++)
	{
		int k2=0;
		while(k2<a)
		{
			Pb[k1*Nsite+k2]=-Pb[k1*Nsite+k2]; //sign changes
			k2++;
		}
	}

	for(int k2=0; k2<Nsite; k2++)	//add the second column
	{
		if(k2==j)
			Pa[Nsite*(Nfermion+1)+k2]=1.0;
		else
			Pa[Nsite*(Nfermion+1)+k2]=0.0;
		if(k2==a)
			Pb[Nsite*(Nfermion+1)+k2]=1.0;
		else
			Pb[Nsite*(Nfermion+1)+k2]=0.0;
	}

	//Calling lapack's dgemm() to do the matrix product
	ldc=m=n=Nfermion+2;
	lda=ldb=k=Nsite;
	alpha=1.0;
	beta=0.0;

	zgemm("C", "N", &m, &n, &k, &alpha, Pa, &lda, Pb, &ldb, &beta, C, &ldc);

	zgetrf(&m, &n, C, &ldc, piv, &info);

	if(info!=0)
		rho=0.0;
	else
		for(int h=0; h<(Nfermion+2);h++)
		{
			if((piv[h]-1)!=h)	
				rho=-rho*C[h*(Nfermion+2)+h];
			else
				rho=rho*C[h*(Nfermion+2)+h];
		}

	free(Pa);
	free(Pb);
	free(C);
	free(piv);

	return rho;
}



// this calculates the two point green's function
//
void TwoPntGreenMx(cdouble * g_2) //green_2 is Nsite*Nsite in dimension
{
	int flag=0;

	for(int i=0; i<Nsite; i++) {
		for(int j=0; j<Nsite; j++) {
			g_2[i*Nsite+j]=0.0;
		}
	}

	for(int i=0; i<Nsite; i++) {
		for(int j=i; j<Nsite; j++) {
//			if(flag%isize==irank) {	//when the job is parallelized, this is part of MPI to recongize each CPU.
			g_2[i*Nsite+j]=TwoPntGreen(i, j);
			if(i!=j) {
				g_2[j*Nsite+i]=conj(g_2[i*Nsite+j]);
			}
		} // end for(int j=i; j<Nsitel j++)
//			} // end if(flag%isize=irank) if we want to include gr
// Parallel processing
//
/*		#ifdef MYMPI

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
		}  // end for(int j=i; j<Nsitel j++), if we want to include parallel processing


#ifdef MYMPI
	MPI::COMM_WORLD.Bcast(g_2, Nsite*Nsite, MPI::DOUBLE_COMPLEX, 0);
#endif
//	cout<<"rank"<<irank<<": the flag of rho matrix is "<<flag<<endl;

} // end TwoPntGreenMatrix with parallel processing

*/
/************************************************************En*d*******************************************************************/
