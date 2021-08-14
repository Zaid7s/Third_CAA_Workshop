F90      = gfortran
F90FLAGS = -g -c -std=f2008 -O -W -fcheck=all -Wline-truncation -Wuninitialized -fbounds-check -Wconversion -Wsurprising -Wpedantic -Warray-temporaries -Wconversion-extra -Wextra -ffpe-trap='zero','overflow','underflow','invalid' -fbacktrace -static
LDFLAGS  = 
LIBS     =

COMPILE_OBJ = $(F90) $(F90FLAGS)

LINK_AND_LOAD = $(F90) $(LDFLAGS) $(LIBS)

CMD = run.exe

OBJS = Main.o \
	   GetMeanValues.o \
	   GetMean_And_Pet.o \
	   GetAreaDerivative.o \
	   GetJacobians.o \
	   GetLUSGS.o \
	   GetDX.o \
	   GetiMaxValue.o \
       Compressible1DFunctions.o \
	   GetExactSol.o \
	   GetLMatrix.o \
	   GetUMatrix.o \
	   GetDMatrix.o \
	   GetRHSVector.o \
	   GetDissipation.o \
	   GetdEdXVector.o \
	   GetInflowBC.o \
	   GetInflowBC_Mean_And_Pet.o \
	   GetDeltaA.o \
	   GetSpatialNOdQdT.o \
	   GetSpatialDerivative.o \
	   GetOutflowBC.o \
	   GetCFL_Ramping.o \
	   GetCFL_Ramping2.o \
	   GetOutflowBC_Mean_And_Pet.o \
	   ForwardSweep.o \
	   BackwardSweep.o \
	   GetRungeKutta.o \
	   GetTriDi.o \
	   GetDiag.o \
	   fft1d.o \
	   Precision_Def.o \
	   BlockTriDiSolver1D.o

$(CMD): $(OBJS)
	($(LINK_AND_LOAD) -o $(@) $(OBJS));
	(./compile.sh)

Main.o: \
 Main.f90 \
 GetMeanValues.o \
 GetMean_And_Pet.o	
	$(COMPILE_OBJ) Main.f90

GetMeanValues.o: \
 GetMeanValues.f90 \
 GetLUSGS.o	\
 GetJacobians.o	\
 GetTriDi.o	\
 GetDX.o	\
 GetCFL_Ramping.o	\
 GetCFL_Ramping2.o	\
 GetiMaxValue.o	\
 GetExactSol.o	\
 GetAreaDerivative.o	
	$(COMPILE_OBJ) GetMeanValues.f90

GetMean_And_Pet.o: \
 GetMean_And_Pet.f90 \
 GetLUSGS.o	\
 GetJacobians.o	\
 GetTriDi.o	\
 GetDX.o	\
 GetiMaxValue.o	\
 GetExactSol.o	\
 fft1d.o \
 Precision_Def.o \
 GetAreaDerivative.o	
	$(COMPILE_OBJ) GetMean_And_Pet.f90

GetiMaxValue.o: \
 GetiMaxValue.f90
	$(COMPILE_OBJ) GetiMaxValue.f90

GetDX.o: \
 GetDX.f90
	$(COMPILE_OBJ) GetDX.f90

GetExactSol.o: \
 GetExactSol.f90 \
 Compressible1DFunctions.o
	$(COMPILE_OBJ) GetExactSol.f90

Compressible1DFunctions.o: \
 Compressible1DFunctions.f90
	$(COMPILE_OBJ) Compressible1DFunctions.f90

GetAreaDerivative.o: \
 GetAreaDerivative.f90
	$(COMPILE_OBJ) GetAreaDerivative.f90

GetJacobians.o: \
 GetJacobians.f90
	$(COMPILE_OBJ) GetJacobians.f90

GetLUSGS.o: \
 GetLUSGS.f90 \
 GetRHSVector.o \
 GetRungeKutta.o \
 GetLMatrix.o \
 GetDMatrix.o \
 GetUMatrix.o \
 ForwardSweep.o \
 BackwardSweep.o 
	$(COMPILE_OBJ) GetLUSGS.f90

GetLMatrix.o: \
 GetLMatrix.f90
	$(COMPILE_OBJ) GetLMatrix.f90

GetUMatrix.o: \
 GetUMatrix.f90
	$(COMPILE_OBJ) GetUMatrix.f90

GetDMatrix.o: \
 GetDMatrix.f90
	$(COMPILE_OBJ) GetDMatrix.f90

GetRHSVector.o: \
 GetRHSVector.f90 \
 GetdEdXVector.o \
 GetSpatialNOdQdT.o \
 GetDissipation.o
	$(COMPILE_OBJ) GetRHSVector.f90

GetSpatialNOdQdT.o: \
 GetSpatialNOdQdT.f90
	$(COMPILE_OBJ) GetSpatialNOdQdT.f90

GetDissipation.o: \
 GetDissipation.f90
	$(COMPILE_OBJ) GetDissipation.f90

GetdEdXVector.o: \
 GetdEdXVector.f90 \
 GetSpatialDerivative.o	\
 GetInflowBC.o \
 GetInflowBC_Mean_And_Pet.o \
 GetOutflowBC.o \
 GetOutflowBC_Mean_And_Pet.o
	$(COMPILE_OBJ) GetdEdXVector.f90

GetInflowBC.o: \
 GetInflowBC.f90 \
 GetDeltaA.o
	$(COMPILE_OBJ) GetInflowBC.f90

GetDeltaA.o: \
 GetDeltaA.f90
	$(COMPILE_OBJ) GetDeltaA.f90

GetInflowBC_Mean_And_Pet.o: \
 GetInflowBC_Mean_And_Pet.f90 \
 GetDeltaA.o
	$(COMPILE_OBJ) GetInflowBC_Mean_And_Pet.f90

GetSpatialDerivative.o: \
 GetSpatialDerivative.f90
	$(COMPILE_OBJ) GetSpatialDerivative.f90

GetOutflowBC.o: \
 GetOutflowBC.f90
	$(COMPILE_OBJ) GetOutflowBC.f90

GetOutflowBC_Mean_And_Pet.o: \
 GetOutflowBC_Mean_And_Pet.f90
	$(COMPILE_OBJ) GetOutflowBC_Mean_And_Pet.f90

ForwardSweep.o: \
 ForwardSweep.f90
	$(COMPILE_OBJ) ForwardSweep.f90

BackwardSweep.o: \
 BackwardSweep.f90
	$(COMPILE_OBJ) BackwardSweep.f90

GetRungeKutta.o: \
 GetRungeKutta.f90 \
 GetdEdXVector.o \
 GetSpatialNOdQdT.o \
 GetDissipation.o
	$(COMPILE_OBJ) GetRungeKutta.f90

fft1d.o:   \
 fft1d.f90 \
 Precision_Def.o
	$(COMPILE_OBJ) fft1d.f90

Precision_Def.o: \
 Precision_Def.f90
	$(COMPILE_OBJ)  Precision_Def.f90

GetCFL_Ramping.o: \
 GetCFL_Ramping.f90
	$(COMPILE_OBJ) GetCFL_Ramping.f90

GetCFL_Ramping2.o: \
 GetCFL_Ramping2.f90
	$(COMPILE_OBJ) GetCFL_Ramping2.f90

GetTriDi.o: \
 GetTriDi.f90 \
 GetRHSVector.o	\
 GetRungeKutta.o \
 GetDiag.o \
 BlockTriDiSolver1D.o
	$(COMPILE_OBJ) GetTriDi.f90

GetDiag.o: \
 GetDiag.f90
	$(COMPILE_OBJ) GetDiag.f90

BlockTriDiSolver1D.o: \
 BlockTriDiSolver1D.f90
	$(COMPILE_OBJ) BlockTriDiSolver1D.f90

clean:
	rm *.o; rm *.mod
