######################  HEAD DEF  #################################
CC= mpicxx
MPIFLAG= -DMYMPI
LIBS= -L/usr/usc/mkl/10.2.1/lib/em64t/ -lmkl -lmkl_lapack64
#-i-dynamic -I/usr/usc/mkl/10.2.1/include \
      #	-lmkl -lmkl_lapack64
OBJ= NoiseCorrelationTrap.o Correlations_Linux.o Greens_Linux.o
TARGET= noiseAmp2.out

$(TARGET): $(OBJ)
	$(CC) -o $(TARGET) $(OBJ) $(MPIFLAG) $(LIBS)

%.o: %.cpp stdafx.h
	$(CC) -c $< $(MPIFLAG) $(LIBS)

.PHONY: clean
clean:
	rm *.o
#######################  END HEAD  #################################
