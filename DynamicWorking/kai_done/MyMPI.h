#ifdef MYMPI    //input data sharing through MPI

#include "stdafxDynKatie.h" //this is my header based on his
//it has all the variables mentioned below, but it doesn't define them
//or say anything about them. So right now this is kinda pointless.

	MPI::Status istatus; //we found this before main()
	//in the noisecorrelation etc.
        
	//we found all of this in the middle of main()
        MPI::COMM_WORLD.Bcast(&Nfermion, 1, MPI::INT, 0);
        MPI::COMM_WORLD.Bcast(&Nsite, 1, MPI::INT, 0);
        MPI::COMM_WORLD.Bcast(&T, 1, MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(&V, 1, MPI::DOUBLE, 0);
	MPI::COMM_WORLD.Bcast(&U, 1, MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(&ratioA, 1, MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(&ratioB, 1, MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(&PeriodicFlag, 1, MPI::CHAR, 0);
        MPI::COMM_WORLD.Bcast(&Nbzone, 1, MPI::INT, 0);
        MPI::COMM_WORLD.Bcast(&BorS, 1, MPI::CHAR, 0);
        MPI::COMM_WORLD.Bcast(&K, 1, MPI::INT, 0);

	MPI::COMM_WORLD.Bcast(&xdisp, 1, MPI::DOUBLE, 0);
	MPI::COMM_WORLD.Bcast(&itime, 1, MPI::INT, 0);
	MPI::COMM_WORLD.Bcast(&ftime, 1, MPI::INT, 0);
	MPI::COMM_WORLD.Bcast(&tstep, 1, MPI::INT, 0);

#endif
